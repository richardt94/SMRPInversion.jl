cd(@__DIR__)
using SMRPInversion
using HiQGA.transD_GP, PyPlot
##
L = 50. # length of side of square in m
ϕ = 13π/36 # magnetic field inclination in radians.
θ = 0 # horizontal angle between loop orientation and magnetic north, in radians.
# pulse moments used for the experiment, in Ampere-seconds.
qgrid = [0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] 
Be = 50000*1e-9 # Earth field magnitude in Tesla.
σ = [0.001] # Earth conductivity assumed
t = Vector{Float64}() # thickness of each layer, empty if only a halfspace

## water sauration discretization
zstart = 1.
extendfrac, dz = 1.028, 1
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=50, showplot=true)
## define water sat model
w = 0.01 * ones(length(zboundaries))
w[(zboundaries .> 30) .& (zboundaries .< 45)] .= 0.4
##
linearsat = true # linear saturation
amponly = true # amplitude only, or invert for phase
mult = false # additive or multiplicative noise
noise_mle = false # maximum likelihood estimate for noise
sounding = SMRPInversion.create_synthetic(w, σ, t, Be, ϕ, L, zboundaries, qgrid; noise_mle, amponly,
                noise_frac=0.05, mult, linearsat)
##
@info "@info χ² is $(2*transD_GP.get_misfit(w, sounding))"