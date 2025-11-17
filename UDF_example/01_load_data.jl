using HiQGA, PyPlot, MAT, DelimitedFiles
cd(@__DIR__)
using SMRPInversion
## Load the processed data from a GMR FID sounding
example_fid = matread("/home/547/ar0754/z67/SMR/UDF_SMR_inversion_ready/SMR_UD01/FID_70ms.mat")
t = example_fid["time_fid"][:]
fid_qt = example_fid["coil_1_fid"]
##
V0, ϕ = SMRPInversion.get_sounding_curve(t, fid_qt)
## params for the sounding - some of these are stored in the MATLAB file
# but others (e.g. field inclination) are from site info for the survey
# or separate ASCII files
q = example_fid["pulse_moment"]
freq = example_fid["detect_frequency"]
L = 50 # 100/2 I think, half side of square length
inclination = -63.35 * π/180 # degrees to radians
θ = 9.13 *π/180  # declination in radians
Be = 2π*freq/SMRPInversion.γh
resist_data = readdlm("/home/547/ar0754/z67/SMR/UDF_SMR_inversion_ready/SMR_UD01/SMR_UD01_res_profile.txt")
c = 1 ./ resist_data[2:end,1] # conductivity
thick = Vector{Float64}(resist_data[2:end-1,2]) # conductivity model thicknesses
σt = SMRPInversion.ConductivityModel(c,thick)
## define a depth grid for the modelling and inversion
zstart = 0.5
extendfrac, dz = 1.028, 1.
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=70, showplot=true)
##
linearsat = true
amponly = true
mult = false
phaserev = false
##
phaserev && (ϕ=-ϕ)
F = SMRPInversion.MRSForward_square(L, zboundaries, q[:], inclination, θ, Be, σt)
sounding = SMRPInversion.newSMRSounding(V0[:], ϕ[:], F, linearsat=linearsat, amponly=amponly, mult=mult, showplot=true)
##
