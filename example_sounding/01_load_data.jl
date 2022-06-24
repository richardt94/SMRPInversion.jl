using transD_GP, PyPlot, MAT, DelimitedFiles, Revise
cd(@__DIR__)
includet("../SMRPI.jl")
includet("../ProcessingTools.jl")
## Load the processed data from a GMR FID sounding
example_fid = matread("../../../example_data/FID_40ms.mat")
t = example_fid["time_fid"][:]
fid_qt = example_fid["coil_1_fid"]
##
V0, ϕ = ProcessingTools.get_sounding_curve(t, fid_qt)
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
inclination = 43.9 * π/180 #degrees to radians
θ = 0 #loop oriented mag. north
Be = 2π*freq/SMRPI.γh
resist_data = readdlm("../../../example_data/example_res_profile.txt")
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
phaserev = true
##
phaserev && (ϕ=-ϕ)
F = SMRPI.MRSForward_square(L, zboundaries, q[:], inclination, 0, Be, σt)
sounding = SMRPI.newSMRSounding(V0[:], ϕ[:], F, linearsat=linearsat, amponly=amponly, mult=mult, showplot=true)
##
GMR_res = readdlm("../../../example_data/conductive_earth_inversion_FID_40ms/conductive_earth_inversion_FID_40ms_1d_inversion.txt")
wc_gmr = GMR_res[1:end-1,4]
z_gmr = (GMR_res[1:end-1, 1] .+ GMR_res[1:end-1, 1])/2
