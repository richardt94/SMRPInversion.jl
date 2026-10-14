# HierarchicalNoise.jl — Paper 1 centrepiece infrastructure.
#
# Adds a hierarchical heteroscedastic noise model to the SMR inversion:
#
#     σᵢ² = s · σ̂ᵢ²
#
# where σ̂ᵢ² is a per-gate empirical variance shape (from stack scatter on
# field data, or a synthetic recipe), and s is a scalar hyperparameter sampled
# jointly with the water-content model. The log-normal prior on s has
# per-run kwargs `s_prior_center` and `s_prior_logwidth` so the E3 money
# figure can vary the prior centre to demonstrate posterior invariance.
#
# The centrepiece contribution vs Taylor 2022 SMRPI:
#   • Taylor's SMRSoundingUnknown profiles a SINGLE scalar σ via ML per sample
#     → biased estimator, no shape, no prior.
#   • This module keeps σ̂ᵢ² SHAPE fixed (heteroscedastic) and SAMPLES the
#     scale s jointly → unbiased posterior on s, invariant to prior centring
#     when the data are informative.
#
# HiQGA wiring:
#   • s enters HiQGA as a nuisance parameter (log s, uniform bounds).
#   • Convention: amponly → nuvec = [log s].
#                 complex → nuvec = [offset_ϕ, log10(stretch_ϕ), log s].
#   • The log-normal prior on s is imposed by ADDING -log P(s) to get_misfit,
#     because HiQGA nuisance sampling assumes uniform priors within bounds.
#
# Likelihood (amponly branch; the term that costs 80 % of the correctness risk
# is the log-det one):
#
#     -log L(w, s | d) = ½ Σᵢ rᵢ² / (s·σ̂ᵢ²) + ½ Σᵢ log(s·σ̂ᵢ²) + const
#                     = ½ Σᵢ rᵢ² / (s·σ̂ᵢ²) + (N/2) log s + ½ Σ log σ̂ᵢ² + const
#
#   plus the log-normal prior contribution
#
#     -log P(s) = ½ (log s − μ)² / τ² + const  ,  μ = log(s_prior_center),
#                                                τ = s_prior_logwidth.
#
# The Σ log σ̂ᵢ² term is constant in (w, s) and dropped from the get_misfit
# return.
#
# Validation gates (must all pass before any E-rerun uses this module):
#   V1: fix w to truth, sample s only → posterior brackets truth, width ≈ analytic.
#   V2: two independent T=1 chains → s posteriors agree.
#   V3: prior-centre sweep ×0.5/×2 on one T1 case → s posteriors converge to
#       same place.
# See `test/test_hierarchical_noise.jl`.

using LinearAlgebra
import HiQGA.transD_GP.AbstractOperator.get_misfit

export SMRSoundingHierarchical, newSMRSoundingHierarchical,
       make_synthetic_hierarchical

"""
    SMRSoundingHierarchical <: SMRSounding

Sounding with heteroscedastic noise model σᵢ² = s · σ̂ᵢ², where σ̂ᵢ² is a fixed
per-gate empirical variance shape and s is sampled jointly with the model.

Fields:
- `V0`, `ϕ`, `Fm`   — as `SMRSoundingKnown`.
- `σ̂²_V0`, `σ̂²_ϕ`  — per-gate empirical variance shapes (amplitude, phase).
                       `σ̂²_ϕ = Float64[]` when `amponly = true`.
- `s`               — current scale (mutable; set by `extractnu!`).
- `s_prior_μ`       — log-normal prior mean, `= log(s_prior_center)`.
- `s_prior_τ`       — log-normal prior stddev.
- other fields identical to `SMRSoundingKnown`.
"""
mutable struct SMRSoundingHierarchical <: SMRSounding
    V0        :: Vector{<:Real}
    ϕ         :: Vector{<:Real}
    Fm        :: SNMRForward.MRSForward
    σ̂²_V0     :: Vector{<:Real}
    σ̂²_ϕ      :: Vector{<:Real}
    s         :: Float64
    s_prior_μ :: Float64
    s_prior_τ :: Float64
    linearsat :: Bool
    amponly   :: Bool
    offset_ϕ  :: Real
    stretch_ϕ :: Real
end

"""
    newSMRSoundingHierarchical(V0, ϕ, Fm;
                                σ̂_V0, σ̂_ϕ = nothing,
                                s_prior_center = 1.0,
                                s_prior_logwidth = 1.0,
                                s_init = 1.0,
                                linearsat = true, amponly = true,
                                offset_ϕ = 0.0, stretch_ϕ = 1.0)

Convenience constructor. Squares σ̂_V0 (and σ̂_ϕ if provided) into the internal
variance fields. `s_prior_center` and `s_prior_logwidth` set the log-normal
prior on s — the two knobs the E3 money figure varies.

`s_init` is the starting value used by the sampler before nuisance updates;
does not affect the posterior at convergence.
"""
function newSMRSoundingHierarchical(V0, ϕ, Fm;
                                    σ̂_V0::AbstractVector,
                                    σ̂_ϕ = nothing,
                                    s_prior_center::Real = 1.0,
                                    s_prior_logwidth::Real = 1.0,
                                    s_init::Real = 1.0,
                                    linearsat::Bool = true,
                                    amponly::Bool = true,
                                    offset_ϕ::Real = 0.0,
                                    stretch_ϕ::Real = 1.0)
    length(V0) == length(ϕ) ||
        throw(ArgumentError("V0 and ϕ must have same length"))
    length(σ̂_V0) == length(V0) ||
        throw(ArgumentError("σ̂_V0 must have same length as V0"))
    amponly && (isapprox(offset_ϕ, 0.0) && isapprox(stretch_ϕ, 1.0)) ||
        amponly && throw(ArgumentError("amponly requires no phase modification"))
    if !amponly
        (σ̂_ϕ === nothing) &&
            throw(ArgumentError("complex-inversion hierarchical requires σ̂_ϕ"))
        length(σ̂_ϕ) == length(ϕ) ||
            throw(ArgumentError("σ̂_ϕ must have same length as ϕ"))
    end
    s_prior_center > 0 || throw(ArgumentError("s_prior_center must be > 0"))
    s_prior_logwidth > 0 || throw(ArgumentError("s_prior_logwidth must be > 0"))
    s_init > 0 || throw(ArgumentError("s_init must be > 0"))

    σ̂²_V0 = Float64.(σ̂_V0) .^ 2
    σ̂²_ϕ  = amponly ? Float64[] : Float64.(σ̂_ϕ) .^ 2

    SMRSoundingHierarchical(V0, ϕ, Fm, σ̂²_V0, σ̂²_ϕ,
                            Float64(s_init),
                            log(Float64(s_prior_center)),
                            Float64(s_prior_logwidth),
                            linearsat, amponly, offset_ϕ, stretch_ϕ)
end

# ── HiQGA nuisance extraction ───────────────────────────────────────────────

"""
    extractnu!(S::SMRSoundingHierarchical, nuvec)

Nuisance vector convention:
  • amponly = true:  nuvec = [log(s)]
  • amponly = false: nuvec = [offset_ϕ, log10(stretch_ϕ), log(s)]

HiQGA samples the nuisance vector uniformly within its bounds; bounds on
log(s) should be wide (e.g. [log(0.01), log(100)]) so the log-normal prior on
s does the shaping.
"""
function extractnu!(S::SMRSoundingHierarchical, nuvec::AbstractVector)
    if S.amponly
        length(nuvec) == 1 ||
            throw(ArgumentError("amp-only hierarchical: nuvec must be length 1 [log s]"))
        S.s = exp(nuvec[1])
    else
        length(nuvec) == 3 ||
            throw(ArgumentError("complex hierarchical: nuvec must be length 3 [offset_ϕ, log10(stretch_ϕ), log s]"))
        S.offset_ϕ = nuvec[1]
        S.stretch_ϕ = 10.0 ^ nuvec[2]
        S.s = exp(nuvec[3])
    end
    nothing
end

# ── Likelihood ──────────────────────────────────────────────────────────────

"""
    get_misfit(w::Vector, S::SMRSoundingHierarchical)

Returns -log P(w, s | d) up to a w,s-independent constant:

    ½ Σ rᵢ²/(s σ̂ᵢ²)  +  (N/2) log s  +  ½ (log s − μ)²/τ²

The N/2 log s term is what prevents s inflating freely; profiling optimises
it away, sampling must keep it.
"""
function get_misfit(w::Vector{<:Real}, S::SMRSoundingHierarchical)
    Vres, ϕres = applynu(S, w)
    s = S.s
    s > 0 || return Inf

    res_amp = S.V0 .- Vres
    chi2 = sum(res_amp .^ 2 ./ (s .* S.σ̂²_V0))
    Ndata = length(res_amp)

    if !S.amponly
        res_ϕ = rem2pi.(S.ϕ .- ϕres, RoundNearest)
        chi2 += sum(res_ϕ .^ 2 ./ (s .* S.σ̂²_ϕ))
        Ndata += length(res_ϕ)
    end

    log_det = Ndata * log(s)                        # +Σ log σ̂ᵢ² dropped (const)
    log_s = log(s)
    prior = (log_s - S.s_prior_μ)^2 / (S.s_prior_τ ^ 2)

    return 0.5 * (chi2 + log_det + prior)
end

# ── Convenience: build a synthetic with true noise + hier prior ─────────────

"""
    make_synthetic_hierarchical(w, σ, t, Be, ϕ_inc, R, zgrid, qgrid;
                                 noise_additive_nV = 50.0,
                                 σ̂_V0_nV = nothing,
                                 s_prior_center = 1.0,
                                 s_prior_logwidth = 1.0,
                                 s_init = 1.0,
                                 rseed = 131,
                                 mult = false, linearsat = true, amponly = true,
                                 square = false, θ = 0.0)

Draw a noisy sounding with additive Gaussian noise of `noise_additive_nV` per
gate (uniform shape by default), then wrap it in an `SMRSoundingHierarchical`
where the sampler's σ̂ᵢ² SHAPE matches the true noise shape.

`σ̂_V0_nV` — if `nothing`, defaults to a uniform vector at `noise_additive_nV`
(the truth). If a vector is passed, that is the shape (in nV) presented to the
sampler; use this to construct wrong-shape studies for E4 panel c-hier.

Returns `(S::SMRSoundingHierarchical, true_noise_V::Real, V0_obs::Vector)`.
"""
function make_synthetic_hierarchical(w, σ, t, Be, ϕ_inc, R, zgrid, qgrid;
                                     noise_additive_nV::Real = 50.0,
                                     σ̂_V0_nV = nothing,
                                     s_prior_center::Real = 1.0,
                                     s_prior_logwidth::Real = 1.0,
                                     s_init::Real = 1.0,
                                     rseed::Int = 131,
                                     mult::Bool = false,
                                     linearsat::Bool = true,
                                     amponly::Bool = true,
                                     square::Bool = false,
                                     θ::Real = 0.0)
    noise_V = noise_additive_nV * 1e-9
    S_known = create_synthetic(w, σ, t, Be, ϕ_inc, R, zgrid, qgrid;
                               noise_additive = noise_V, mult = mult,
                               amponly = amponly, linearsat = linearsat,
                               noise_mle = false, showplot = false,
                               rseed = rseed, square = square, θ = θ)
    Nq = length(S_known.V0)
    σ̂_V0 = (σ̂_V0_nV === nothing) ? fill(noise_V, Nq) : (σ̂_V0_nV .* 1e-9)
    S = newSMRSoundingHierarchical(
        S_known.V0, S_known.ϕ, S_known.Fm;
        σ̂_V0 = σ̂_V0, σ̂_ϕ = amponly ? nothing : (σ̂_V0 ./ max.(S_known.V0, 1e-12)),
        s_prior_center = s_prior_center, s_prior_logwidth = s_prior_logwidth,
        s_init = s_init,
        linearsat = linearsat, amponly = amponly)
    (S = S, true_noise_V = noise_V, V0_obs = S_known.V0)
end
