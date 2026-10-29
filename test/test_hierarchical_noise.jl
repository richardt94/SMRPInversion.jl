# test_hierarchical_noise.jl — validation gates V1 / V2 / V3 for the paper's
# hierarchical heteroscedastic noise model (σᵢ² = s · σ̂ᵢ²) added in
# `src/HierarchicalNoise.jl`.
#
# These tests must all pass before any E-rerun uses the hierarchical model.
# V1: fixed simple model, only-s sampled → posterior brackets truth, width ≈ analytic.
# V2: two independent T=1 chains → s posteriors agree.
# V3: prior-centre sweep ×0.5/×2 → s posteriors converge to same place.

using Test, Statistics
using SMRPInversion, SNMRForward
using HiQGA: transD_GP

# ────────────────────────────────────────────────────────────────────────────
# Analytic reference for V1
# ────────────────────────────────────────────────────────────────────────────
#
# Fix w = w_true. Only s is sampled. Data is d = f(w_true) + noise with noise
# ~ N(0, σ̂² I). Since w is fixed, residuals r = d - f(w_true) = noise, so
# Σ r²/σ̂² ~ χ²_N.
#
# The get_misfit returns
#   M(s) = ½ Σ r²/(s σ̂²) + (N/2) log s + ½ (log s − μ)² / τ²
#
# Setting dM/ds = 0 gives the log posterior mode:
#   0 = -½ Σ r²/(s² σ̂²) + N/(2s) + (log s − μ) / (s τ²)
#   s_mode ≈ (Σ r²/σ̂²) / N   when the prior is broad (τ → ∞)
#
# So with the true noise realisation R² = Σ r²/σ̂², the sampler should recover
# s_mode ≈ R²/N, and R²/N → 1 as N grows (χ²_N/N → 1).
#
# We use N_data = 24 gates → expected posterior s to sit near 1.0, with
# spread ~ sqrt(2/N_data) ≈ 0.29 in relative units under a broad prior.

# ── Common setup (small, fast case) ─────────────────────────────────────────
function build_test_case(; rseed = 42)
    zall, znall, zbnd = transD_GP.setupz(1.0, 1.028; dz = 1.0, n = 40, showplot = false)
    w_true = fill(0.02, length(zbnd))
    w_true[(zbnd .>= 20.0) .& (zbnd .< 30.0)] .= 0.25
    σ_layers = [0.05, 0.05]
    t_layers = [500.0]
    qgrid = 10.0 .^ range(log10(0.1), log10(15.0), length = 24) |> collect
    Be    = 50_000e-9
    ϕ_inc = deg2rad(-60.0)
    R     = 50.0
    (w_true = w_true, zbnd = zbnd, σ_layers = σ_layers, t_layers = t_layers,
     qgrid = qgrid, Be = Be, ϕ_inc = ϕ_inc, R = R, rseed = rseed)
end

# ────────────────────────────────────────────────────────────────────────────
@testset "hierarchical sounding construction + extractnu!" begin
    tc = build_test_case()
    r = SMRPInversion.make_synthetic_hierarchical(
        tc.w_true, tc.σ_layers, tc.t_layers, tc.Be, tc.ϕ_inc, tc.R,
        tc.zbnd, tc.qgrid;
        noise_additive_nV = 50.0,
        s_prior_center = 1.0, s_prior_logwidth = 1.0,
        s_init = 1.0, rseed = tc.rseed)

    S = r.S
    @test isa(S, SMRSoundingHierarchical)
    @test S.amponly == true
    @test length(S.σ̂²_V0) == length(tc.qgrid)
    @test all(S.σ̂²_V0 .> 0)
    @test S.s == 1.0

    # extractnu! for amponly = 1-element nuisance
    SMRPInversion.extractnu!(S, [log(2.5)])
    @test S.s ≈ 2.5
    SMRPInversion.extractnu!(S, [log(0.4)])
    @test S.s ≈ 0.4

    # Wrong-length nuvec should throw
    @test_throws ArgumentError SMRPInversion.extractnu!(S, [0.0, 0.0])
end

# ────────────────────────────────────────────────────────────────────────────
@testset "get_misfit likelihood terms — analytic sanity" begin
    tc = build_test_case(rseed = 7)
    r = SMRPInversion.make_synthetic_hierarchical(
        tc.w_true, tc.σ_layers, tc.t_layers, tc.Be, tc.ϕ_inc, tc.R,
        tc.zbnd, tc.qgrid;
        noise_additive_nV = 50.0,
        s_prior_center = 1.0, s_prior_logwidth = 1.0e6,   # essentially flat prior
        rseed = tc.rseed)
    S = r.S

    # Compute residuals directly at w_true
    Vres, _ = SMRPInversion.applynu(S, tc.w_true)
    res = S.V0 .- Vres
    N = length(res)

    # ∂M/∂(log s) at s: derivative w.r.t. log s of
    #   M(s) = ½ Σ r²/(s σ̂²) + (N/2) log s + tiny prior
    # is:  -½/s · Σ r²/σ̂² + N/2 + (log s − μ)/τ²
    # Under flat prior (τ→∞): setting to zero gives
    #   s_star = Σ r²/σ̂² / N
    R² = sum(res.^2 ./ S.σ̂²_V0)
    s_star = R² / N

    # M at s_star should be a local minimum. Check that derivative at s = s_star
    # in log-s space is small (< 1e-4).
    Δ = 1e-6
    SMRPInversion.extractnu!(S, [log(s_star) - Δ])
    M_minus = SMRPInversion.get_misfit(tc.w_true, S)
    SMRPInversion.extractnu!(S, [log(s_star) + Δ])
    M_plus  = SMRPInversion.get_misfit(tc.w_true, S)
    dM_dlogs = (M_plus - M_minus) / (2Δ)
    @test isapprox(dM_dlogs, 0.0; atol = 1e-3)

    # Curvature at the mode should be positive (minimum, not saddle)
    SMRPInversion.extractnu!(S, [log(s_star)])
    M0 = SMRPInversion.get_misfit(tc.w_true, S)
    @test M_minus > M0
    @test M_plus > M0

    # Reasonable magnitude: since w is truth and noise realisation is fair,
    # s_star should be within a factor 2 of 1.0 for N=24.
    @test 0.5 < s_star < 2.0
end

# ────────────────────────────────────────────────────────────────────────────
@testset "get_misfit prior term shape" begin
    tc = build_test_case()
    r = SMRPInversion.make_synthetic_hierarchical(
        tc.w_true, tc.σ_layers, tc.t_layers, tc.Be, tc.ϕ_inc, tc.R,
        tc.zbnd, tc.qgrid;
        noise_additive_nV = 50.0,
        s_prior_center = 1.0, s_prior_logwidth = 0.5,
        rseed = 5)
    S = r.S

    # With a narrow log-normal prior (τ = 0.5), M(s) should curve up sharply
    # away from s = 1.0 even holding likelihood fixed. Test that misfit at
    # s = e^2 is much larger than at s = 1.
    SMRPInversion.extractnu!(S, [0.0])
    M_at_1 = SMRPInversion.get_misfit(tc.w_true, S)
    SMRPInversion.extractnu!(S, [2.0])
    M_at_e2 = SMRPInversion.get_misfit(tc.w_true, S)
    @test M_at_e2 - M_at_1 > 5.0        # prior term ≈ (2)²/(2·0.25) = 8
end

# ────────────────────────────────────────────────────────────────────────────
@testset "wrong-length nuisance / bad kwargs throw" begin
    tc = build_test_case()
    r = SMRPInversion.make_synthetic_hierarchical(
        tc.w_true, tc.σ_layers, tc.t_layers, tc.Be, tc.ϕ_inc, tc.R,
        tc.zbnd, tc.qgrid;
        rseed = tc.rseed)
    S = r.S

    @test_throws ArgumentError SMRPInversion.extractnu!(S, [0.0, 0.0, 0.0])   # length 3 in amponly
    @test_throws ArgumentError SMRPInversion.newSMRSoundingHierarchical(
        S.V0, S.ϕ, S.Fm; σ̂_V0 = [1.0, 2.0])                                  # wrong length
    @test_throws ArgumentError SMRPInversion.newSMRSoundingHierarchical(
        S.V0, S.ϕ, S.Fm; σ̂_V0 = fill(1e-9, length(S.V0)), s_prior_center = -1.0)
    @test_throws ArgumentError SMRPInversion.newSMRSoundingHierarchical(
        S.V0, S.ϕ, S.Fm; σ̂_V0 = fill(1e-9, length(S.V0)), s_prior_logwidth = 0.0)
end

# NOTE: V1, V2, V3 full MCMC-based validation gates require running trans_D
# with OptionsNuisance for a hierarchical sounding. This is orchestrated from
# the paper1/ synthetics repo (see `paper1/hierarchical_validation.jl`) so it
# can use the 4-worker setup and PT ladder. Unit tests above cover the maths
# and interface; the sampler-integration tests live alongside the E-series
# scripts.
