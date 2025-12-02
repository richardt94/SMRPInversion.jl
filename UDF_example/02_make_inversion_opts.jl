## make options for the stationary GP
fileprefix = "SMR_example_"*(linearsat ? "linearsat_" : "logsat_")*
                            (amponly ? "amponly_" : "phase_")*
                            (mult ? "multnoise_" : "")
K = transD_GP.GP.OrstUhn()
nmin, nmax = 2, 40 # max number of basis functions 
linearbounds = [.001 0.35] # saturation bounds
fbounds = linearsat ? linearbounds : log10.(linearbounds)
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])] # perturbation for depth
sdev_prop = 0.07*diff(fbounds, dims=2)[:] # perturbation for property
sdev_dc = 0.008*diff(fbounds, dims=2)[:] # perturbation for DC shifts in property
xall = permutedims(collect(1:1.0:length(zboundaries))) # depths modelled at in number of layers
xbounds = permutedims([extrema(xall)...]) # depth bounds in number of layers
λ, δ = [2], (linearsat ? fbounds[1]/10 : 1e-3) # correlation length, and nugget in saturation
## Initialize a stationary GP using these options
using Random
Random.seed!(12)
opt = transD_GP.OptionsStat(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        fdataname = fileprefix,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        save_freq = 50,
                        sdev_dc = sdev_dc,
                        quasimultid = false,
                        K = K
                        )
