# examples/occam_on_udf_data.jl
#
# End-to-end example: load a UDF-format GMR sounding, extract the sounding
# curve, build a `SMRSoundingKnown`, and run the deterministic Occam
# inversion (`SMRPInversion.deterministic_invert`).
#
# The .mat FID file and the associated resistivity profile live on NCI
# (`/g/data/z67/...`) and are not part of this repo. If the paths below
# don't exist locally, this script prints a clear message and exits — so
# it can safely sit alongside working examples without breaking anything.

using SMRPInversion, SNMRForward
using HiQGA: transD_GP
using MAT, DelimitedFiles, PyPlot, Statistics, Printf

# Adjust these two paths to point at your UDF-format data:
const MAT_PATH = "/g/data/z67/ar0754/z67/SMR/UDF_SMR_inversion_ready/SMR_UD01/FID_70ms.mat"
const RES_PATH = "/g/data/z67/ar0754/z67/SMR/UDF_SMR_inversion_ready/SMR_UD01/SMR_UD01_res_profile.txt"

if !(isfile(MAT_PATH) && isfile(RES_PATH))
    @warn "UDF example data not available at $MAT_PATH — skipping."
    exit()
end

# ── Load .mat and extract sounding curve ────────────────────────────────────
d      = matread(MAT_PATH)
t      = vec(d["time_fid"])
fid_qt = d["coil_1_fid"]
q_raw  = vec(d["pulse_moment"])
freq   = float(d["detect_frequency"][1])

V0_est, ϕ_est = SMRPInversion.get_sounding_curve(t, fid_qt)
V0 = vec(V0_est); ϕ = vec(ϕ_est)

# Truncate the first few pulse moments to avoid dead-time contamination
ntrunc = 4
q  = q_raw[1+ntrunc:end]
V0 = V0[1+ntrunc:end]
ϕ  = ϕ[1+ntrunc:end]

# ── Load ERT-derived conductivity model ─────────────────────────────────────
resist_data = readdlm(RES_PATH)
c     = 1.0 ./ resist_data[2:end, 1]
thick = Float64.(resist_data[2:end-1, 2])
σt    = SMRPInversion.ConductivityModel(c, thick)

# ── Site geometry (UD01, per meta_UDF.txt) ─────────────────────────────────
# 100 m square loop → L = 50 (half-side per SNMRForward convention)
L      = 50.0
inc    = deg2rad(-63.35)
θ      = deg2rad(9.13)
Be     = 2π * freq / SMRPInversion.γh

# ── Depth grid (matches SMRPI conventions) ─────────────────────────────────
zstart, extendfrac, dz = 1.0, 1.028, 1.0
zall, znall, zbnd = transD_GP.setupz(zstart, extendfrac; dz = dz, n = 70,
                                     showplot = false)

# ── Build the sounding with a fixed σ_V₀ so Occam has a known noise level ──
# For real data, pick σ_V₀ from a per-site derivation (fit residual,
# late-time envelope, or MCMC's noise-MLE posterior). Here we hard-code
# 30 nV as a placeholder — replace with a defensible per-site number.
σ_V = 30e-9
F = SMRPInversion.MRSForward_square(L, zbnd, q, inc, θ, Be, σt)
sounding = SMRPInversion.newSMRSounding(V0, ϕ, F;
    σ_V0 = fill(σ_V, length(V0)),
    σ_ϕ  = σ_V ./ V0,
    linearsat = true, amponly = true, mult = false, showplot = false)

# ── Occam inversion ────────────────────────────────────────────────────────
op, mnew, χ², λ², oidx = SMRPInversion.deterministic_invert(sounding;
    regtype = :R1)
w_occam = mnew[end][oidx[end]]
χ²_final = χ²[end][oidx[end]]
λ²_final = λ²[end][oidx[end]][1]

@printf("UDF UD01 Occam inversion complete.\n")
@printf("  χ² at convergence: %.2f  (target %d)\n", χ²_final, length(V0))
@printf("  log10(λ²) chosen : %.2f\n", log10(λ²_final))
@printf("  peak WC          : %.3f at z = %.1f m\n",
        maximum(w_occam), zbnd[argmax(w_occam)])

# ── Plot recovered profile ─────────────────────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(6, 7))
ax.step(w_occam, zbnd, color="darkorange", linewidth=1.7,
        label="Occam R1 (χ²=$(round(χ²_final, digits=2)))")
ax.set_ylim(120, 0); ax.set_xlim(0, 0.35)
ax.set_xlabel("Water content"); ax.set_ylabel("Depth (mbgl)")
ax.grid(true, alpha=0.3); ax.legend(fontsize=9, loc="lower right")
ax.set_title("UDF UD01 — deterministic Occam recovery", fontsize=10)
fig.tight_layout()
savefig("examples/occam_on_udf_UD01.png", dpi=140, bbox_inches="tight")
plt.close("all")
@info "wrote examples/occam_on_udf_UD01.png"
