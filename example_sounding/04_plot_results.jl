using Statistics
##
close("all")
transD_GP.getchi2forall(opt)
##
istothepow = false
@assert !(linearsat & istothepow)
opt.xall[:] .= zboundaries
transD_GP.plot_posterior(sounding, opt, burninfrac=0.5, figsize=(10,6), qp1=0.05, qp2=0.95, nbins=50, istothepow=istothepow, cmappdf="inferno", CIcolor=["c", "b"], fsize=12,
    vmaxpc=1, pdfnormalize=true, plotmean=false, lwidth=1)
ax = gcf().axes
linearsat ? ax[1].set_xlabel("fractional water content") : ax[1].set_xlabel("log\$_{10}\$ water content") 
linearsat || ax[1].plot(istothepow ? wc_gmr : log10.(wc_gmr), z_gmr, "w-")
linearsat && ax[1].plot(wc_gmr, z_gmr, "w-")
savefig(fileprefix*"post.png", dpi=600)
## nuisance histograms
if !amponly 
    transD_GP.plot_posterior(sounding, optn, burninfrac=0.5, nbins=50, figsize=(5,4))
    savefig(fileprefix*"nuisance.png", dpi=600)
end    
## swarm plots
if amponly
    SMRPI.plot_model_field(sounding, opt, decfactor=10, lcolor="k", modelalpha=0.08)
else    
    SMRPI.plot_model_field(sounding, opt, optn, decfactor=10, lcolor="k", modelalpha=0.08)
end
savefig(fileprefix*"swarm.png", dpi=600)
## noise estimates
ndata = amponly ? length(sounding.V0) : 2*length(sounding.V0)
F = transD_GP.assembleTat1(opt, :U, temperaturenum=1)
est_σ2 = exp.(2/ndata * F)/ndata
est_σ = sqrt.(est_σ2)
@info mean(est_σ)
