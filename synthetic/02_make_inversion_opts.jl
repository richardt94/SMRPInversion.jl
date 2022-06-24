## make options for the stationary GP
fileprefix = "SMR_example_"*(linearsat ? "linearsat_" : "logsat_")*
                            (amponly ? "amponly_" : "")*
                            (mult ? "multnoise_" : "")*
                            (noise_mle ? "mle_" : "")
K = transD_GP.GP.OrstUhn()
nmin, nmax = 2, 40
linearbounds = [.001 0.5]
fbounds = linearsat ? linearbounds : log10.(linearbounds)
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.07*diff(fbounds, dims=2)[:]
demean = false
sdev_dc = 0.008*diff(fbounds, dims=2)[:]
sampledc = true
xall = permutedims(collect(1:1.0:length(zboundaries)))
xbounds = permutedims([extrema(xall)...])
λ, δ = [2], (linearsat ? 0.01 : 0.0001)
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
                        demean = demean,
                        sampledc = sampledc,
                        sdev_dc = sdev_dc,
                        quasimultid = false,
                        K = K
                        )
