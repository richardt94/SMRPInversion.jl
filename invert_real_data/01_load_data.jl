using HiQGA, PyPlot, MAT, DelimitedFiles, Revise
cd(@__DIR__)
# includet("../SMRPI.jl")
using SMRPInversion
# includet("../ProcessingTools.jl")
## Load the processed data from a GMR FID sounding
example_fid = matread("/g/data/z67/rlt118/SNMR_data/OK18_38/FID_80ms_noemi.mat")
t = example_fid["time_fid"][:]
fid_qt = example_fid["coil_1_fid"]
##
V0, ϕ = SMRPInversion.get_sounding_curve(t, fid_qt)
## params for the sounding - some of these are stored in the MATLAB file
# but others (e.g. field inclination) are from site info for the survey
# or separate ASCII files
q = example_fid["pulse_moment"]
freq = example_fid["detect_frequency"]
## truncate dataset
ntrunc = 4
q = q[:,1+ntrunc:end]
V0 = V0[:,1+ntrunc:end]
ϕ = ϕ[:, 1+ntrunc:end]
##
L = 50
inclination = 44.2 * π/180 #degrees to radians
θ = 0 #loop oriented mag. north
Be = 2π*freq/SMRPI.γh
resist_data = readdlm("/g/data/z67/rlt118/SNMR_data/OK18_38/OK18_38_res_profile.txt")
c = 1 ./ resist_data[2:end,1]
thick = Vector{Float64}(resist_data[2:end-1,2])
σt = SMRPI.ConductivityModel(c,thick)
## define a depth grid for the modelling and inversion
zstart = 1.
extendfrac, dz = 1.028, 1
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=70, showplot=true)
gcf()
##
linearsat = true
amponly = true
mult = false
phaserev = false
##
phaserev && (ϕ=-ϕ)
F = SMRPI.MRSForward_square(L, zboundaries, q[:], inclination, 0, Be, σt)
sounding = SMRPI.newSMRSounding(V0[:], ϕ[:], F, linearsat=linearsat, amponly=amponly, mult=mult, showplot=true)
##
