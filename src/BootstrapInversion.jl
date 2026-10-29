# Bootstrap uncertainty for Occam-style deterministic SMR inversion.
#
# Two flavours:
#   :parametric  — resample fresh Gaussian noise from the asserted σ_V0 and
#                  add to the noise-free forward at the best-fit w0. Standard
#                  parametric bootstrap; matches the noise model deterministic
#                  Occam and MCMC already assume.
#   :residual    — resample the observed residuals `V0 - |K̃·w0|` with
#                  replacement. Non-parametric; useful when the noise
#                  distribution is uncertain and you have field data.
#
# Both produce an ensemble of Occam solutions. Their per-depth spread is a
# non-linear, non-Gaussian uncertainty band directly comparable to
# `linearised_posterior_covariance` (the Gaussian ±2σ approximation) and to
# an MCMC p5-p95 credible interval. Where the three agree, the linearised
# assumption is safe; where bootstrap is wider, non-linearity or multi-modality
# matters and the linearised band is over-confident.
#
# Cost: one Occam per replicate. Occam converges in a few seconds per SMR
# sounding, so N = 100-200 replicates typically fits in a few minutes.

using LinearAlgebra, SparseArrays, Random, Statistics

export bootstrap_invert

"""
    bootstrap_invert(S::SMRSoundingKnown;
                     N = 200,
                     w0 = nothing,
                     method = :parametric,
                     rseed = nothing,
                     verbose = false,
                     kwargs...) -> NamedTuple

Bootstrap uncertainty for `deterministic_invert`.

Arguments:
- `S`      : sounding with known noise (`SMRSoundingKnown`).
- `N`      : number of bootstrap replicates. Default 200.
- `w0`     : centre model. Required for `:parametric` (noise is added to the
             forward at `w0`) and for `:residual` (residuals are computed
             against the forward at `w0`). If `nothing` (default), first
             solves the unperturbed problem with `deterministic_invert` and
             uses that solution. Pass a pre-computed best fit to save one
             inversion.
- `method` : `:parametric` — add fresh Gaussian noise ∼ N(0, σ_V0²) to
             |K̃·w0| per replicate. Matches the noise model MCMC / linearised
             assume, so a fair three-way comparison.
             `:residual` — draw the observed residuals `V0 - |K̃·w0|` with
             replacement (non-parametric). Useful for field data where the
             noise is not obviously Gaussian.
             Default `:parametric`.
- `rseed`  : base seed for reproducibility. If given, replicate `i` uses
             `MersenneTwister(rseed + i)`. If `nothing`, uses the global RNG.
- `verbose`: log per-replicate progress. Default `false`.
- `kwargs...` : forwarded to `deterministic_invert` (`wstart`, `wref`,
                `regtype`, `β²`, `lo`, `hi`, `λ²min`, `λ²max`, `target`, etc.).
                **Critical**: pass the SAME `wstart` / `wref` / `regtype` /
                `β²` you used for the original inversion. Bootstrap
                characterises noise sensitivity around a fixed regularisation
                choice, not sensitivity to the prior itself. To probe wref
                sensitivity, use a B3-style reference-model sweep instead.

Returns a NamedTuple with:
- `w0`      : centre model used (either supplied or from the initial fit).
- `w_reps`  : `(nz × N)` matrix of replicate Occam solutions. Columns
              corresponding to failed inversions are all `NaN`.
- `n_ok`    : number of replicates that converged.
- `q05, q50, q95` : per-cell quantiles across successful replicates.
- `mean_w, std_w` : per-cell mean and standard deviation across replicates.
- `λ²_reps` : final `λ²` chosen by each replicate's line search — spread
              here is the regularisation-weight sensitivity the linearised
              approach implicitly misses.
- `χ²_reps` : final `χ²` achieved by each replicate.
- `method`  : the bootstrap method used, echoed.

Notes:
- Bootstrap does not capture `wref` uncertainty (each replicate uses the
  same reference model). Combine with a B3-style `wref` sweep if you need
  that.
- For amplitude-only inversion (the only mode `SMRDetOperator` supports),
  the noise-free forward is `|K̃·w0|` and bootstrap perturbations are added
  to that amplitude. Very small perturbations that push individual gates
  negative are clipped to a tiny positive floor before inversion, as
  negative amplitudes are unphysical.
- Each replicate builds its own `SMRDetOperator` (a fresh trapezoidal-
  weighted kernel copy). Marginal cost — dominated by the inversion itself.
"""
function bootstrap_invert(S::SMRSoundingKnown;
        N::Int = 200,
        w0::Union{Nothing, AbstractVector} = nothing,
        method::Symbol = :parametric,
        rseed::Union{Nothing, Int} = nothing,
        verbose::Bool = false,
        kwargs...)

    method in (:parametric, :residual) ||
        throw(ArgumentError("method must be :parametric or :residual, got $method"))

    # Solve the unperturbed problem if the caller didn't supply a centre model
    if w0 === nothing
        verbose && @info "bootstrap_invert: solving unperturbed Occam for w0"
        op0, mnew0, _, _, oidx0 = deterministic_invert(S; kwargs...)
        w0 = mnew0[end][oidx0[end]]
    end
    w0 = collect(Float64, w0)

    # Compute the noise-free forward at w0 via the amplitude operator.
    # Use S.amponly to confirm the amplitude assumption matches the operator.
    op_ref = SMRDetOperator(S)
    length(w0) == size(op_ref.K̃, 2) ||
        throw(ArgumentError("w0 length ($(length(w0))) must equal kernel column count ($(size(op_ref.K̃, 2)))"))
    V_hat = abs.(op_ref.K̃ * w0)

    # Residuals only needed for :residual method
    residuals = method === :residual ? (S.V0 .- V_hat) : Float64[]

    nz = length(w0)
    w_reps  = fill(NaN, nz, N)
    λ²_reps = fill(NaN, N)
    χ²_reps = fill(NaN, N)
    n_ok = 0

    for i in 1:N
        rng = rseed === nothing ? Random.GLOBAL_RNG : MersenneTwister(rseed + i)

        # Generate a bootstrap V0
        if method === :parametric
            noise_i = randn(rng, length(V_hat)) .* S.σ_V0
        else
            idx = rand(rng, 1:length(residuals), length(residuals))
            noise_i = residuals[idx]
        end
        V_boot = max.(V_hat .+ noise_i, 1e-30)  # clip: negative amplitudes are unphysical

        # Build the bootstrap sounding — swap only V0, keep everything else
        S_boot = SMRSoundingKnown(V_boot, S.ϕ, S.σ_V0, S.σ_ϕ, S.Fm,
                                   S.linearsat, S.amponly, S.offset_ϕ, S.stretch_ϕ)

        try
            _, mnew_i, χ²hist_i, λ²hist_i, oidx_i = deterministic_invert(S_boot; kwargs...)
            w_reps[:, i] = mnew_i[end][oidx_i[end]]
            λ²_reps[i]   = λ²hist_i[end][oidx_i[end]][1]
            χ²_reps[i]   = χ²hist_i[end][oidx_i[end]]
            n_ok += 1
            verbose && (i % 25 == 0) && @info "bootstrap_invert: replicate $i / $N done ($n_ok OK)"
        catch e
            verbose && @warn "bootstrap replicate $i failed: $e"
        end
    end

    # Per-cell quantiles + mean/std, robust to NaN rows
    q05    = _rowwise(w_reps, x -> quantile(x, 0.05))
    q50    = _rowwise(w_reps, x -> quantile(x, 0.50))
    q95    = _rowwise(w_reps, x -> quantile(x, 0.95))
    mean_w = _rowwise(w_reps, mean)
    std_w  = _rowwise(w_reps, std)

    return (w0 = w0,
            w_reps  = w_reps,
            n_ok    = n_ok,
            q05     = q05, q50 = q50, q95 = q95,
            mean_w  = mean_w, std_w = std_w,
            λ²_reps = λ²_reps, χ²_reps = χ²_reps,
            method  = method)
end

# NaN-robust per-row reduction: apply `f` to the non-NaN entries of each row,
# returning NaN for all-NaN rows. Used for the summary stats above.
function _rowwise(mat::AbstractMatrix{Float64}, f)
    nz = size(mat, 1)
    out = Vector{Float64}(undef, nz)
    @inbounds for i in 1:nz
        row = filter(!isnan, view(mat, i, :))
        out[i] = isempty(row) ? NaN : f(row)
    end
    out
end
