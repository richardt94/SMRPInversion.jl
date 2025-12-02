## set up McMC
using Distributed
nsamples, nchains, nchainsatone = 10_001, 4, 1
Tmax = 2.5
addprocs(nchains)
@info "workers are $(workers())"
@everywhere using Distributed, HiQGA.transD_GP, SMRPInversion
## run McMC
@time begin
        transD_GP.main(opt, sounding; Tmax, nsamples, nchains, nchainsatone)
end        
## close the worker pool
rmprocs(workers())