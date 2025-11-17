using Statistics
##
close("all")
transD_GP.getchi2forall(opt)
##
opt.xall[:] .= zboundaries
transD_GP.plot_posterior(sounding, opt, burninfrac=0.5, qp1=0.05, qp2=0.95, nbins=50,
    vmaxpc=1.0, pdfnormalize=false, plotmean=false, lwidth=1)
ax = gcf().axes
ax[1].set_xlabel("fractional water content")
ax[1].step(zboundaries, "w-")
## amplitude swarm plots
if sounding.amponly
    SMRPInversion.plot_model_field(sounding, opt, decfactor=10, lcolor="k", modelalpha=0.08)
else    
    SMRPInversion.plot_model_field(sounding, opt, optn, decfactor=10, lcolor="k", modelalpha=0.08)
end  
## noise estimates
if noise_mle
    ndata = amponly ? length(sounding.V0) : 2*length(sounding.V0)
    F = transD_GP.assembleTat1(opt, :U, temperaturenum=1)
    est_σ2 = exp.(2/ndata * F)/ndata
    est_σ = sqrt.(est_σ2)
    @info mean(est_σ)
end