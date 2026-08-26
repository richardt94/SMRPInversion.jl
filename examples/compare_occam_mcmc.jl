# examples/compare_occam_mcmc.jl
#
# Side-by-side comparison of the deterministic Occam inversion
# (SMRPInversion.deterministic_invert) and the probabilistic
# trans-dimensional MCMC (HiQGA.transD_GP.main) on the same synthetic
# amplitude-only SMR sounding.
#
# Not part of `test/runtests.jl` — MCMC is expensive and needs distributed
# workers, both of which are awkward in a CI context. Run manually:
#
#     julia --project=. examples/compare_occam_mcmc.jl
#
# The MCMC settings here are deliberately modest (100 001 samples × 4
# chains) so this finishes in a minute or two on a laptop. For a paper
# figure you'd want more samples.

using Distributed
const NCHAINS = 4
addprocs(max(0, NCHAINS - nworkers()))
@everywhere using Distributed, HiQGA.transD_GP, SMRPInversion

using SNMRForward, PyPlot, Statistics

# ── Synthetic aquifer at 40–50 m, 25 % WC, resistive halfspace ─────────────
zall, znall, zbnd = transD_GP.setupz(1.0, 1.028; dz=1.0, n=60, showplot=false)
w_true = fill(0.02, length(zbnd))
w_true[(zbnd .>= 40.0) .& (zbnd .< 50.0)] .= 0.25

qgrid  = 10.0 .^ range(log10(0.1), log10(15.0), length=24) |> collect
L      = 50.0    # 100 m loop → half-side
φ_inc  = deg2rad(-60.0)
Be     = 50_000e-9
σ_lay  = [0.05, 0.05]; t_lay = [500.0]

sounding = SMRPInversion.create_synthetic(
    w_true, σ_lay, t_lay, Be, φ_inc, L, zbnd, qgrid;
    noise_additive = 50e-9, mult = false, amponly = true, linearsat = true,
    noise_mle = false, showplot = false, rseed = 42)

@info "Synthetic sounding: peak V0 = $(round(1e9*maximum(sounding.V0), digits=1)) nV"

# ── Occam ──────────────────────────────────────────────────────────────────
n_z    = length(zbnd)
wref   = fill(0.02, n_z)
wstart = fill(0.02, n_z)
op, mnew, χ², λ², oidx = SMRPInversion.deterministic_invert(
    sounding; wstart = wstart, wref = wref, regtype = :R1)
w_occam = mnew[end][oidx[end]]
χ²_occam = χ²[end][oidx[end]]
@info "Occam: χ² = $(round(χ²_occam, digits=2)), peak w = $(round(maximum(w_occam), digits=3))"

# ── MCMC (via HiQGA trans-D) ───────────────────────────────────────────────
K            = transD_GP.GP.OrstUhn()
fbounds      = [0.001 0.5]
sdev_pos     = [0.05 * abs(diff([extrema(znall)...])[1])]
sdev_prop    = 0.07 * diff(fbounds, dims=2)[:]
sdev_dc      = 0.008 * diff(fbounds, dims=2)[:]
xbounds      = [1.0  Float64(n_z)]
xall_fresh   = permutedims(collect(1:1.0:n_z))

opt = transD_GP.OptionsStat(
    nmin = 2, nmax = 40,
    xbounds = xbounds, fbounds = fbounds,
    fdataname = "compare_occam_mcmc_",
    xall = xall_fresh,
    λ = [4.0], δ = 0.01,
    sdev_prop = sdev_prop, sdev_pos = sdev_pos,
    save_freq = 50, sdev_dc = sdev_dc,
    quasimultid = false, K = K,
)
transD_GP.main(opt, sounding; Tmax = 2.5,
    nsamples = 100_001, nchains = NCHAINS, nchainsatone = 1)
rmprocs(workers())

# Extract MCMC posterior quantiles
opt.xall[:] .= zbnd
F = transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
samples = hcat([f[:] for f in F]...)
q05 = mapslices(x -> quantile(x, 0.05), samples; dims = 2)[:]
q50 = mapslices(x -> quantile(x, 0.50), samples; dims = 2)[:]
q95 = mapslices(x -> quantile(x, 0.95), samples; dims = 2)[:]
@info "MCMC: peak posterior median = $(round(maximum(q50), digits=3))"

# ── Comparison figure ──────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(6, 7))
ax.fill_betweenx(zbnd, q05, q95, alpha=0.30, color="steelblue", label="MCMC p5–p95")
ax.step(q50, zbnd, color="steelblue", linewidth=1.5, label="MCMC median")
ax.step(w_occam, zbnd, color="darkorange", linewidth=1.7,
        label="Occam R1 (χ²=$(round(χ²_occam, digits=1)))")
ax.step(w_true, zbnd, "k-", linewidth=1.5, label="True")
ax.set_ylim(120, 0); ax.set_xlim(0, 0.35)
ax.set_xlabel("Water content"); ax.set_ylabel("Depth (mbgl)")
ax.grid(true, alpha=0.3); ax.legend(fontsize=9, loc="lower right")
ax.set_title("Occam vs MCMC on a synthetic 10 m aquifer at 40 m", fontsize=10)
fig.tight_layout()
savefig("examples/compare_occam_mcmc.png", dpi=140, bbox_inches="tight")
plt.close("all")
@info "wrote examples/compare_occam_mcmc.png"
