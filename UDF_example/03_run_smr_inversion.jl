## set up McMC
using Distributed
nsamples, nchains, nchainsatone = 100_001, 4, 1
Tmax = 2.5
addprocs(nchains)
@info "workers are $(workers())"
@everywhere using Distributed, HiQGA.transD_GP, SMRPInversion
## run McMC
@time begin
    if sounding.amponly
        transD_GP.main(opt, sounding, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
    else    
        transD_GP.main(opt, optn, sounding, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
    end
end        
## close the worker pool
rmprocs(workers())