module ProcessingTools

using FFTW, LinearAlgebra, Statistics

export get_sounding_curve

function analytic_signal(signal)
    sigf = fft(signal)
    freq = fftfreq(length(signal))
    anaf = ((freq .>= 0) .+ (freq .> 0)) .* sigf
    ifft(anaf)
end

#the estimate functions do provide error estimates but
# they do not account for uncertainties in the
# preprocessing (filtering etc.) - be careful
# using those estimates
function est_E0_T2(t, timeseries)
    amp = abs.(analytic_signal(timeseries))
    log_amp = log.(amp)

    G = hcat(-t, ones(length(t)))
    regr = pinv(G)
    est = regr*log_amp
    T2 = 1/est[1]
    E0 = exp(est[2])
    modelled_amp = E0 * exp.(-t/T2)
    vard = var(amp .- modelled_amp)
    σlog = sqrt(vard) ./ amp
    Cd = diagm(σlog.^2)
    Cm = regr * Cd * regr'
    σm = sqrt.(diag(Cm))
    σT2 = T2^2 * σm[1]
    σE0 = E0*σm[2]

    (E0, T2, σE0, σT2)
end

function unwrap!(phase)
    for i in 2:length(phase)
        while phase[i] - phase[i-1] >= pi
            phase[i] -= 2pi
        end
        while phase[i] - phase[i-1] <= -pi
            phase[i] += 2pi
        end
    end
end

function est_ω_ϕ(t, timeseries)
    arg = angle.(analytic_signal(timeseries))
    unwrap!(arg)
    G = hcat(t, ones(length(t)))
    regr = pinv(G)
    ω, ϕ = regr*arg
    vard = var(ω*t .+ ϕ .- arg)
    Cm = vard * regr * regr'
    σω, σϕ = sqrt.(diag(Cm))

    (ω, ϕ, σω, σϕ)
end

function full_estimates(t, ts)
    (est_E0_T2(t,ts)..., est_ω_ϕ(t,ts)...)
end

function get_sounding_curve(t, fid_qt)
    est = mapslices(x -> full_estimates(t, x), fid_qt; dims=1)
    getindex(i) = map(x->x[i], est)

    getindex(1), getindex(6)
end

end