# Deterministic (Occam-style Gauss-Newton) SMR inversion.
#
# Wraps an `SMRSoundingKnown` into a HiQGA `Operator1D` that
# `HiQGA.transD_GP.gradientinv` can consume. The SMR forward is
#     V(q) = K̃(q, z) · w(z)
# — linear in water content `w`, complex-valued in general. This module
# implements the amplitude-only inversion, |V(q)| = |K̃·w|, using the
# analytical amplitude-Jacobian
#
#     J[j, i] = ∂|V_j| / ∂w_i = Re(conj(V_j) · K̃[j, i]) / |V_j|.
#
# The Occam framework (smoothest model that fits the data to within
# χ² = n_data) follows Constable, Parker & Constable (1987), with
# HiQGA's `gradientinv` providing the outer line-search over λ² and the
# Gauss-Newton mechanics.

using LinearAlgebra, SparseArrays

export SMRDetOperator, deterministic_invert, linearised_posterior_covariance

"""
    SMRDetOperator <: HiQGA.transD_GP.AbstractOperator.Operator1D

Deterministic-inversion wrapper around an `SMRSoundingKnown`, satisfying the
interface HiQGA's Occam-style `gradientinv` expects.

Fields:
- `S`      : the underlying SMR sounding (data, noise, forward object).
- `K̃`      : `(nq × nz)` complex kernel matrix such that `V = K̃ * w`
             (SNMRForward kernel × trapezoidal quadrature weights on `zgrid`).
- `ρ`      : latest evaluated model iterate (kept in the field HiQGA expects).
- `nfixed` : number of fixed depth cells at the top (always 0 here).
- `J`      : `(nq × nz)` real amplitude Jacobian, populated in-place by
             `getresidual(op, w; computeJ=true)`.
- `W`      : `diag(1 ./ σ_V0)`, sparse — used for chi-square weighting.
- `res`    : `(nq,)` current residual `|V(w)| − V0`, updated by
             `getresidual`.
"""
mutable struct SMRDetOperator <: transD_GP.AbstractOperator.Operator1D
    S      :: SMRSoundingKnown
    K̃      :: Matrix{ComplexF64}
    ρ      :: Vector{Float64}
    nfixed :: Int
    J      :: Matrix{Float64}
    W      :: SparseMatrixCSC{Float64,Int}
    res    :: Vector{Float64}
end

"""
    SMRDetOperator(S::SMRSoundingKnown) -> SMRDetOperator

Build a deterministic-inversion operator from an amplitude-only SMR sounding.

Precomputes `K̃ = trapezoidal_weights .* kernel'` so the amplitude forward is
`|K̃ * w|`. Throws if `S.amponly` is false — the analytical Jacobian is
amplitude-specific.
"""
function SMRDetOperator(S::SMRSoundingKnown)
    S.amponly || throw(ArgumentError("SMRDetOperator currently supports amplitude-only inversion."))
    kernel = S.Fm.kernel                       # (nz × nq) complex
    dz     = S.Fm.dz                           # (nz-1,) layer thicknesses
    nz, nq = size(kernel)
    nz == length(S.Fm.zgrid) || error("kernel row count must match zgrid length")

    # Trapezoidal quadrature weights at each depth boundary
    wt = zeros(Float64, nz)
    wt[1]   = 0.5 * dz[1]
    wt[end] = 0.5 * dz[end]
    @inbounds for i in 2:nz-1
        wt[i] = 0.5 * (dz[i-1] + dz[i])
    end
    K̃ = permutedims(wt .* kernel)              # (nq × nz), so V = K̃ * w

    J   = zeros(Float64, nq, nz)
    W   = sparse(Diagonal(1.0 ./ S.σ_V0))
    res = zeros(Float64, nq)
    ρ   = zeros(Float64, nz)
    SMRDetOperator(S, K̃, ρ, 0, J, W, res)
end

"""
    getresidual(op::SMRDetOperator, w::Vector{Float64}; computeJ=false)

Evaluate the amplitude-only forward at model `w`, populate `op.res` with the
signed residual `|K̃·w| − V0`, and — if `computeJ=true` — populate `op.J` with
the analytical amplitude Jacobian
`J[j,i] = Re(conj(V[j]) · K̃[j,i]) / |V[j]|`.

Implements `HiQGA.transD_GP.AbstractOperator.getresidual`.
"""
function transD_GP.AbstractOperator.getresidual(op::SMRDetOperator, w::Vector{Float64}; computeJ=false)
    V   = op.K̃ * w
    amp = abs.(V)
    @inbounds for j in eachindex(amp)
        op.res[j] = amp[j] - op.S.V0[j]
    end
    if computeJ
        # Floor on |V_j| prevents NaN when the complex kernel row nearly
        # cancels for a particular w. In that regime the amplitude
        # derivative is genuinely very large (small perturbations to w
        # produce a large relative change in |V|), so a large-but-finite
        # J is the right physical answer.
        @inbounds for i in axes(op.K̃, 2), j in axes(op.K̃, 1)
            denom = max(amp[j], 1e-30)
            op.J[j,i] = real(conj(V[j]) * op.K̃[j,i]) / denom
        end
    end
    op.ρ .= w
    nothing
end

"""
    makereg(r::Symbol, F::SMRDetOperator) -> SparseMatrixCSC

Depth-scaled smoothness operator for `SMRDetOperator`, dispatched over
HiQGA's generic `makereg`.

The base HiQGA `makereg` builds a bare first-difference matrix `Rw =
vcat(0, diff(w))`. On the log-growing depth grid produced by
`transD_GP.setupz(1.0, 1.028; ...)` cells at 100 m are ~15 m thick vs
~1 m at the top, so an unweighted difference penalty over-regularises
the shallow model and under-regularises the deep model. Replaces the
row `i ≥ 2` entries with `±1/√dz[i-1]`, giving `‖Rw‖² ≈ ∫(∂w/∂z)² dz`
— the Constable-Parker-Constable (1987) discretisation-consistent form.

Row 1 is `1/√dz[1]` on the diagonal, anchoring the top cell against
zero at the same 1/√dz scale as the differences. The base `:R1` leaves
row 1 all zeros, so with `β² = 0` (the default) the top cell drifts
freely — problematic for SMR where dead-time truncation makes the top
few metres the least-constrained. For a tight anchor to a non-zero
reference model, set `β² > 0` in `deterministic_invert` and provide a
plausible `wref[1]`.
"""
function transD_GP.makereg(r::Symbol, F::SMRDetOperator)
    dz = F.S.Fm.dz
    nz = size(F.K̃, 2)
    length(dz) == nz - 1 || error("dz length ($(length(dz))) must be nz-1 ($(nz-1))")
    r === :R0 && return sparse(1.0I, nz, nz)
    if r === :R1
        rows = Int[1]; cols = Int[1]; vals = Float64[1.0 / sqrt(dz[1])]
        @inbounds for i in 2:nz
            s = 1.0 / sqrt(dz[i-1])
            push!(rows, i); push!(cols, i-1); push!(vals, -s)
            push!(rows, i); push!(cols, i);   push!(vals,  s)
        end
        return sparse(rows, cols, vals, nz, nz)
    end
    if r === :R2
        R1 = transD_GP.makereg(:R1, F)
        return R1 * R1
    end
    error("unknown regularization: $r")
end

"""
    deterministic_invert(S::SMRSoundingKnown; wstart, wref, kwargs...)
        -> (op, mnew, χ², λ², oidx)

Run Occam-style Gauss-Newton inversion on an amplitude-only SMR sounding.

Uses the analytical Jacobian above; delegates the outer line-search over the
smoothness weight `λ²` to `HiQGA.transD_GP.gradientinv`. Follows Constable,
Parker & Constable (1987) — "smoothest model that fits within χ² target."

Arguments:
- `S`      : sounding with known noise (`SMRSoundingKnown`).
- `wstart` : starting water-content vector, one entry per boundary in
             `S.Fm.zgrid`. Default `0.05` everywhere.
- `wref`   : reference model used only when `β² > 0`. Default `0.05`.
- `lo`, `hi` : per-cell bounds on water content. Default `[1e-4, 0.4]`.
              `hi = 0.4` reflects a plausible upper bound on porosity for
              most sands/silts; tighten further for site-specific priors.
- `regtype` : `:R0` (smallest), `:R1` (first-difference), or `:R2`
              (curvature). `:R1` and `:R2` are dispatched through a
              depth-scaled `makereg` local to this module — see the
              method docstring for `makereg(::Symbol, ::SMRDetOperator)`.
              Default `:R1`.
- `λ²min`, `λ²max` : log10 bounds on the smoothness-weight line search.
                      Default `[-1, 6]`.
- `β²`     : reference-model penalty weight. Default `0` (wref inactive).
             Set `> 0` if you need a tight anchor of the top cell to
             `wref[1]` — the scaled R matrix only weakly anchors it to
             zero on its own.
- `knownvalue` : Gauss-Newton acceptance threshold — a step is accepted
                 if it drops χ² to at most `knownvalue · χ²₀`. Inherited
                 from HiQGA's `gradientinv`. Default `0.7`.
- `nstepsmax` : maximum outer iterations. Default `15`.
- `ntries` : line-search resolution on `λ²`. Default `6`.
- `target` : χ² target. Default `length(S.V0)` (χ²/N = 1).
- Any remaining keyword arguments are forwarded to `gradientinv`.

Returns `(op, mnew, χ², λ², oidx)`:
- `op` is the `SMRDetOperator` — holds the final Jacobian, residual and
  weights. Feed to `linearised_posterior_covariance` for a Gaussian
  uncertainty band around the recovered model.
- The remaining fields are `gradientinv`'s per-iteration histories. The
  final model is `mnew[end][oidx[end]]`.

See also: `linearised_posterior_covariance`.
"""
function deterministic_invert(S::SMRSoundingKnown;
        wstart = fill(0.05, length(S.Fm.zgrid)),
        wref   = fill(0.05, length(S.Fm.zgrid)),
        lo = 1e-4, hi = 0.4,
        regtype = :R1,
        λ²min = -1.0,
        λ²max = 6.0,
        β²    = 0.0,
        knownvalue = 0.7,
        nstepsmax = 15,
        ntries    = 6,
        target    = nothing,
        kwargs...)
    op = SMRDetOperator(S)
    length(wstart) == size(op.K̃, 2) || throw(ArgumentError("wstart length must equal length(S.Fm.zgrid)"))
    length(wref)   == length(wstart) || throw(ArgumentError("wref length must equal wstart length"))
    result = transD_GP.gradientinv(
        collect(Float64, wstart), collect(Float64, wref), op;
        regtype, λ²min, λ²max, β², knownvalue, nstepsmax, ntries, target, lo, hi, kwargs...)
    return (op, result...)
end

"""
    linearised_posterior_covariance(op::SMRDetOperator, λ²::Real;
                                     regtype = :R1)
        -> NamedTuple with fields (Cm, σ, λ²)

Compute the linearised Gaussian posterior covariance of the water-content
model around the point at which `op.J` was last populated.

The formula is the standard linear-inverse result

    C_m = (Jᵀ W² J + λ² · Rᵀ R)⁻¹

where `J` is `op.J`, `W = op.W = diag(1 ./ σ_V₀)` (so `W² = op.W' * op.W` is
the inverse noise covariance), and `R` is the smoothness matrix produced
by `HiQGA.transD_GP.makereg(regtype, op)`. The Gaussian per-depth
uncertainty is `σ = sqrt.(diag(C_m))` — call this on the diagonal after
`deterministic_invert` and you have a ±2σ band directly comparable to
MCMC's p5–p95 credible interval.

**`λ²` is deliberately mandatory.** Omitting the smoothness term
(`λ² = 0`) leaves `Jᵀ W² J`, which is ill-conditioned or singular for the
underdetermined smooth-model inversions this operator is designed for and
produces meaningless uncertainty. Use the final `λ²` chosen by the Occam
line search:

    op, mnew, χ², λ²hist, oidx = deterministic_invert(S; regtype = :R1)
    λ²_final = λ²hist[end][oidx[end]][1]
    result   = linearised_posterior_covariance(op, λ²_final; regtype = :R1)
    σ_cells  = result.σ                     # per-depth 1-σ uncertainty

Keyword `regtype` must match the `regtype` passed to `deterministic_invert`
(default `:R1`).

Notes:
- Assumes `op.J` was populated at the model of interest. `deterministic_invert`
  ends with `op.J` at the final iterate, so most callers don't need to do
  anything explicit. If you re-use the operator for other computations,
  call `HiQGA.transD_GP.AbstractOperator.getresidual(op, w; computeJ = true)`
  first to be safe.
- The returned covariance ignores the [lo, hi] box constraint used during
  inversion. Cells whose mean sits close to a bound have an artificially
  wide Gaussian σ. Clip the ±2σ band to `[lo, hi]` for display if that
  matters.

Returns a NamedTuple with:
- `Cm` — full covariance matrix (dense, `nz × nz`).
- `σ`  — per-cell 1-σ uncertainty, `sqrt.(diag(Cm))`.
- `λ²` — the smoothness weight used (echoed for downstream bookkeeping).
"""
function linearised_posterior_covariance(op::SMRDetOperator, λ²::Real;
                                          regtype::Symbol = :R1)
    R   = transD_GP.makereg(regtype, op)
    JtW = op.J' * op.W
    H   = JtW * JtW' .+ float(λ²) .* Matrix(R' * R)
    Cm  = inv(cholesky(Symmetric(H)))
    σ   = sqrt.(max.(diag(Cm), 0.0))
    return (Cm = Cm, σ = σ, λ² = float(λ²))
end
