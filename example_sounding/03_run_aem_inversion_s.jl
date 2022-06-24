## set up McMC
using Distributed
nsamples, nchains, nchainsatone = 400001, 8, 1 # seems real data needs about 8 chains not 4
Tmax = 2.5
addprocs(nchains)
@info "workers are $(workers())"
@everywhere begin
    using Distributed
    using transD_GP
    !isdefined(@__MODULE__, :SMRPI) && include("../SMRPI.jl")
end 
## run McMC
@time begin
    if sounding.amponly
        transD_GP.main(opt, sounding, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
    else    
        transD_GP.main(opt, optn, sounding, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
    end
end        
## close the worker pool
# rmprocs(workers())
