using HiQGA.transD_GP, PyPlot, Revise
cd(@__DIR__)
includet("../SMRPI.jl")
##
L = sqrt(π*50^2/4)
log_zgrid = 0.5:0.05:2
zgrid = 10 .^ log_zgrid
ϕ = 13π/36
θ = 0
qgrid = [0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
Be = 50000*1e-9#2π*2500/SNMRForward.γh
σ = [0.001]
t = Vector{Float64}()
##
zstart = 1.
extendfrac, dz = 1.028, 1
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=50, showplot=true)
gcf()
##
w = 0.01 * ones(length(zboundaries))
w[(zboundaries .> 30) .& (zboundaries .< 45)] .= 0.4
##
linearsat = true
amponly = false
mult = false
noise_mle = false
sounding = SMRPI.create_synthetic(w, σ, t, Be, ϕ, 50., zboundaries, qgrid, noise_mle=noise_mle, amponly=amponly,
                noise_frac=0.05, mult=mult, linearsat=linearsat)
##
transD_GP.get_misfit(w, sounding)