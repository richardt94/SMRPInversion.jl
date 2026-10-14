# Tests for bootstrap_invert
#
# Covers:
#   1. Parametric bootstrap runs to completion, returns the documented
#      NamedTuple shape, all replicates converge on a well-posed problem.
#   2. Per-cell quantile ordering (q05 <= q50 <= q95) holds.
#   3. λ² spread across replicates is non-degenerate — different noise
#      realisations pick different regularisation strengths (this is the
#      sensitivity that linearised_posterior_covariance omits).
#   4. `rseed` reproducibility.
#   5. Residual-bootstrap variant runs.
#   6. Bootstrap median is close to the centre model w0 (unbiased around
#      the fit).

using Test, Statistics
using SMRPInversion, SNMRForward
using HiQGA: transD_GP

# Reuse the synthetic from test_deterministic.jl
zall, znall, zbnd = transD_GP.setupz(1.0, 1.028; dz=1.0, n=60, showplot=false)
w_true = fill(0.02, length(zbnd))
w_true[(zbnd .>= 20.0) .& (zbnd .< 30.0)] .= 0.25
σ_layers = [0.05, 0.05]
t_layers = [500.0]
qgrid    = 10.0 .^ range(log10(0.1), log10(15.0), length=24) |> collect
sounding = SMRPInversion.create_synthetic(
    w_true, σ_layers, t_layers, 50_000e-9, deg2rad(-60.0), 50.0, zbnd, qgrid;
    noise_additive = 100e-9,
    mult = false, amponly = true, linearsat = true, noise_mle = false,
    showplot = false, rseed = 42,
)

# Best-fit Occam solution for reuse
n = length(zbnd)
wref   = fill(0.02, n)
wstart = fill(0.02, n)
op0, mnew0, χ²0, λ²hist0, oidx0 = SMRPInversion.deterministic_invert(
    sounding; wstart = wstart, wref = wref, regtype = :R1)
w0 = mnew0[end][oidx0[end]]

# ────────────────────────────────────────────────────────────────────────────
@testset "parametric bootstrap shape + convergence" begin
    N = 20
    result = SMRPInversion.bootstrap_invert(
        sounding; N = N, w0 = w0, method = :parametric,
        rseed = 1234,
        wstart = wstart, wref = wref, regtype = :R1)

    # Documented fields present
    @test haskey(result, :w0)
    @test haskey(result, :w_reps)
    @test haskey(result, :n_ok)
    @test haskey(result, :q05) && haskey(result, :q50) && haskey(result, :q95)
    @test haskey(result, :mean_w) && haskey(result, :std_w)
    @test haskey(result, :λ²_reps) && haskey(result, :χ²_reps)
    @test result.method == :parametric

    # Shapes
    @test length(result.w0) == n
    @test size(result.w_reps) == (n, N)
    @test length(result.q05) == n
    @test length(result.λ²_reps) == N

    # All replicates should converge for this well-posed problem
    @test result.n_ok == N
    @test all(!isnan, result.λ²_reps)
    @test all(!isnan, result.χ²_reps)
end

# ────────────────────────────────────────────────────────────────────────────
@testset "quantile ordering" begin
    result = SMRPInversion.bootstrap_invert(
        sounding; N = 30, w0 = w0, method = :parametric,
        rseed = 7, wstart = wstart, wref = wref, regtype = :R1)

    # Per-cell quantile ordering
    for i in 1:n
        @test result.q05[i] <= result.q50[i] <= result.q95[i]
    end

    # Mean should sit inside p5-p95
    for i in 1:n
        @test result.q05[i] <= result.mean_w[i] <= result.q95[i]
    end

    # Standard deviation non-negative and finite
    @test all(result.std_w .>= 0.0)
    @test all(isfinite, result.std_w)
end

# ────────────────────────────────────────────────────────────────────────────
@testset "λ² spread across replicates" begin
    result = SMRPInversion.bootstrap_invert(
        sounding; N = 30, w0 = w0, method = :parametric,
        rseed = 21, wstart = wstart, wref = wref, regtype = :R1)

    # Different noise realisations should pick different λ² — the whole
    # point of bootstrap versus a single-solution linearised covariance.
    λ²s = filter(!isnan, result.λ²_reps)
    @test length(λ²s) >= 20
    @test std(log10.(λ²s)) > 0.1     # at least ~25 % spread in λ²
end

# ────────────────────────────────────────────────────────────────────────────
@testset "rseed reproducibility" begin
    r1 = SMRPInversion.bootstrap_invert(
        sounding; N = 10, w0 = w0, method = :parametric,
        rseed = 999, wstart = wstart, wref = wref, regtype = :R1)
    r2 = SMRPInversion.bootstrap_invert(
        sounding; N = 10, w0 = w0, method = :parametric,
        rseed = 999, wstart = wstart, wref = wref, regtype = :R1)

    @test r1.w_reps == r2.w_reps
    @test r1.λ²_reps == r2.λ²_reps
end

# ────────────────────────────────────────────────────────────────────────────
@testset "residual bootstrap" begin
    result = SMRPInversion.bootstrap_invert(
        sounding; N = 15, w0 = w0, method = :residual,
        rseed = 5, wstart = wstart, wref = wref, regtype = :R1)
    @test result.method == :residual
    @test result.n_ok >= 10    # should mostly converge
    @test all(result.q05 .<= result.q95)
end

# ────────────────────────────────────────────────────────────────────────────
@testset "bootstrap centre near w0" begin
    # With N large enough, the bootstrap MEDIAN should be close to w0 — the
    # bootstrap is characterising noise around the point estimate, not a
    # different estimate.
    result = SMRPInversion.bootstrap_invert(
        sounding; N = 50, w0 = w0, method = :parametric,
        rseed = 314, wstart = wstart, wref = wref, regtype = :R1)

    # Per-cell mean should track w0 within a few percent of the bootstrap
    # standard deviation. Not exact — regularisation strength varies per
    # replicate, biasing the mean slightly toward smoother solutions.
    for i in 1:n
        @test abs(result.mean_w[i] - w0[i]) < max(0.10, 3 * result.std_w[i])
    end
end

# ────────────────────────────────────────────────────────────────────────────
@testset "argument validation" begin
    @test_throws ArgumentError SMRPInversion.bootstrap_invert(
        sounding; N = 2, w0 = w0, method = :bogus,
        wstart = wstart, wref = wref, regtype = :R1)

    # Wrong-length w0
    bad_w0 = fill(0.05, n + 3)
    @test_throws ArgumentError SMRPInversion.bootstrap_invert(
        sounding; N = 2, w0 = bad_w0, method = :parametric,
        wstart = wstart, wref = wref, regtype = :R1)
end

# ────────────────────────────────────────────────────────────────────────────
# Dispatch via deterministic_invert's `uncertainty` kwarg
# ────────────────────────────────────────────────────────────────────────────
@testset "uncertainty kwarg dispatch" begin
    # :none returns nothing in the uncertainty field, everything else intact
    r_none = SMRPInversion.deterministic_invert(sounding;
        wstart = wstart, wref = wref, regtype = :R1, uncertainty = :none)
    @test r_none.uncertainty === nothing
    @test size(r_none.op.J) == (length(sounding.V0), n)   # inversion still ran

    # :laplace returns a Cm + σ + λ² NamedTuple
    r_lap = SMRPInversion.deterministic_invert(sounding;
        wstart = wstart, wref = wref, regtype = :R1, uncertainty = :laplace)
    @test r_lap.uncertainty !== nothing
    @test haskey(r_lap.uncertainty, :Cm) && haskey(r_lap.uncertainty, :σ) && haskey(r_lap.uncertainty, :λ²)
    @test size(r_lap.uncertainty.Cm) == (n, n)
    @test length(r_lap.uncertainty.σ) == n
    @test all(r_lap.uncertainty.σ .>= 0.0)

    # :linearised is an alias for :laplace — should give an equivalent shape
    r_lin = SMRPInversion.deterministic_invert(sounding;
        wstart = wstart, wref = wref, regtype = :R1, uncertainty = :linearised)
    @test r_lin.uncertainty !== nothing
    @test haskey(r_lin.uncertainty, :σ)
    @test r_lin.uncertainty.σ ≈ r_lap.uncertainty.σ

    # :bootstrap returns the bootstrap NamedTuple with quantile arrays
    r_boot = SMRPInversion.deterministic_invert(sounding;
        wstart = wstart, wref = wref, regtype = :R1,
        uncertainty = :bootstrap, N_bootstrap = 10, bootstrap_rseed = 111)
    @test r_boot.uncertainty !== nothing
    @test haskey(r_boot.uncertainty, :w_reps) && haskey(r_boot.uncertainty, :q05)
    @test size(r_boot.uncertainty.w_reps) == (n, 10)
    @test r_boot.uncertainty.method == :parametric
    @test r_boot.uncertainty.n_ok == 10

    # Bogus method throws
    @test_throws ArgumentError SMRPInversion.deterministic_invert(sounding;
        wstart = wstart, wref = wref, regtype = :R1, uncertainty = :bogus)

    # Backward compat: positional unpacking of the first 5 fields still works
    op, mnew, χ², λ², oidx = SMRPInversion.deterministic_invert(sounding;
        wstart = wstart, wref = wref, regtype = :R1, uncertainty = :none)
    @test op isa SMRPInversion.SMRDetOperator
    @test length(mnew) >= 1
end
