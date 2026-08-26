# Tests for SMRDetOperator + deterministic_invert
#
# Covers:
#   1. Analytical amplitude-Jacobian matches a finite-difference reference
#      at a representative model. Numerical validation of the derivation
#      J[j,i] = Re(conj(V_j) · K̃[j,i]) / |V_j|.
#   2. `deterministic_invert` converges on a synthetic aquifer, hitting
#      the χ² = n_data target within reasonable tolerance and recovering
#      the aquifer bump in the right depth range.

using Test, LinearAlgebra, Statistics
using SMRPInversion, SNMRForward
using HiQGA: transD_GP

# ── Reusable synthetic setup ────────────────────────────────────────────────
"""
Build a fresh synthetic sounding with a 10 m aquifer at 20 m depth,
25 % water content, in a resistive half-space. Returns the sounding, the
true water content vector on the inversion grid, and the depth grid.
"""
function synthetic_sounding(; noise_V = 100e-9)
    # depth grid, log-growing to ~200 m
    zall, znall, zbnd = transD_GP.setupz(1.0, 1.028; dz=1.0, n=60, showplot=false)

    # true water content: 25 % between 20 and 30 m, 2 % elsewhere
    w_true = fill(0.02, length(zbnd))
    w_true[(zbnd .>= 20.0) .& (zbnd .< 30.0)] .= 0.25

    # resistive halfspace
    σ_layers = [0.05, 0.05]
    t_layers = [500.0]

    # 24 log-spaced pulse moments 0.1 – 15 A·s
    qgrid = 10.0 .^ range(log10(0.1), log10(15.0), length=24) |> collect

    # 100 m square loop (L is half-side per SNMRForward convention),
    # inclination −60°, Earth field 50 μT
    L, φ_inc, θ_loop, Be = 50.0, deg2rad(-60.0), 0.0, 50_000e-9

    S = SMRPInversion.create_synthetic(
        w_true, σ_layers, t_layers, Be, φ_inc, L, zbnd, qgrid;
        noise_additive = noise_V,
        mult = false, amponly = true, linearsat = true, noise_mle = false,
        showplot = false, rseed = 42,
    )
    return S, w_true, zbnd
end

# Instantiate once for reuse across sub-testsets
sounding, w_true, zbnd = synthetic_sounding()

# ────────────────────────────────────────────────────────────────────────────
@testset "Analytical Jacobian matches finite differences" begin
    op = SMRPInversion.SMRDetOperator(sounding)
    w_eval = copy(w_true)

    # Populate op.J with the analytical Jacobian at w_eval
    transD_GP.AbstractOperator.getresidual(op, w_eval; computeJ = true)
    J_ana = copy(op.J)

    # Central-difference finite Jacobian
    nq, nz = size(J_ana)
    J_fd = zeros(nq, nz)
    h    = 1e-6
    for i in 1:nz
        wp = copy(w_eval); wp[i] += h
        wm = copy(w_eval); wm[i] -= h
        transD_GP.AbstractOperator.getresidual(op, wp; computeJ = false)
        ap = copy(op.res) .+ sounding.V0    # forward |V| = res + V0
        transD_GP.AbstractOperator.getresidual(op, wm; computeJ = false)
        am = copy(op.res) .+ sounding.V0
        J_fd[:, i] .= (ap .- am) ./ (2h)
    end
    # Restore op.J
    transD_GP.AbstractOperator.getresidual(op, w_eval; computeJ = true)

    # Comparison: Frobenius relative error should be < 1e-6 (limited by float64
    # + central-difference step size, not by the analytical formula)
    frob_rel = norm(J_ana .- J_fd) / norm(J_ana)
    @test frob_rel < 1e-6
end

# ────────────────────────────────────────────────────────────────────────────
@testset "deterministic_invert converges to χ² target" begin
    n = length(zbnd)
    wref   = fill(0.02, n)
    wstart = fill(0.02, n)

    op, mnew, χ², λ², oidx = SMRPInversion.deterministic_invert(
        sounding; wstart = wstart, wref = wref, regtype = :R1)

    # final iterate
    wfinal    = mnew[end][oidx[end]]
    χ²_final  = χ²[end][oidx[end]]
    target    = length(sounding.V0)

    # χ² should be within ~10% of target for a well-posed problem
    @test isapprox(χ²_final, target; rtol = 0.1)

    # Solution should live in the physical bounds
    @test all(1e-4 .<= wfinal .<= 0.5)

    # Iterations shouldn't blow up
    @test length(mnew) <= 15
end

# ────────────────────────────────────────────────────────────────────────────
@testset "deterministic_invert recovers aquifer amplitude" begin
    n = length(zbnd)
    wref   = fill(0.02, n)
    wstart = fill(0.02, n)

    op, mnew, χ², λ², oidx = SMRPInversion.deterministic_invert(
        sounding; wstart = wstart, wref = wref, regtype = :R1)
    wfinal = mnew[end][oidx[end]]

    # Aquifer sits at 20–30 m in the truth. Occam's smoothness broadens it,
    # so the recovered peak in a window around the aquifer should be well
    # above the 2 % background but need not match the true 25 % exactly.
    aquifer_mask = (zbnd .>= 15.0) .& (zbnd .<= 40.0)
    peak_in_aquifer = maximum(wfinal[aquifer_mask])
    @test peak_in_aquifer > 0.10
    @test peak_in_aquifer < 0.35

    # Away from the aquifer (below 60 m), we should stay closer to background
    # than to the aquifer amplitude. Occam's smoothness prior lets the aquifer
    # trail into deeper cells, so allow up to ~half the aquifer peak.
    deep_mask = zbnd .> 60.0
    @test maximum(wfinal[deep_mask]) < 0.12
end

# ────────────────────────────────────────────────────────────────────────────
@testset "linearised_posterior_covariance" begin
    n = length(zbnd)
    wref   = fill(0.02, n)
    wstart = fill(0.02, n)
    op, mnew, χ², λ², oidx = SMRPInversion.deterministic_invert(
        sounding; wstart = wstart, wref = wref, regtype = :R1)

    λ²_final = λ²[end][oidx[end]][1]

    # Compute at the same regtype used in inversion
    result = SMRPInversion.linearised_posterior_covariance(op, λ²_final;
                                                            regtype = :R1)

    # Shape sanity
    @test size(result.Cm) == (n, n)
    @test length(result.σ) == n
    @test result.λ² ≈ float(λ²_final)

    # Uncertainty is non-negative and finite everywhere
    @test all(isfinite, result.σ)
    @test all(result.σ .>= 0.0)

    # Kernel sensitivity dies with depth, so σ should be larger deep than
    # shallow (SMR resolvability limit). Use the top/bottom 25 % of the
    # grid so this holds regardless of grid extent.
    n_split = length(zbnd) ÷ 4
    shallow_σ = maximum(result.σ[1:n_split])
    deep_σ    = maximum(result.σ[end - n_split + 1:end])
    @test deep_σ > shallow_σ

    # Increasing λ² tightens the posterior (more smoothing → smaller σ per
    # cell where the smoothness regulariser is active)
    result_looser  = SMRPInversion.linearised_posterior_covariance(
        op, λ²_final / 10; regtype = :R1)
    result_tighter = SMRPInversion.linearised_posterior_covariance(
        op, λ²_final * 10; regtype = :R1)
    @test mean(result_tighter.σ) < mean(result_looser.σ)
end
