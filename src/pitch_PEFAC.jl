#=
Implementation of
Gonzalez, S., & Brookes, M. (2014). PEFAC-a pitch estimation algorithm robust to high levels of noise. IEEE/ACM Transactions on Audio, Speech, and Language Processing, 22(2), 518-530.
=#




pitch_estimation_filter(q::Real, γ::Real, K::Int) =  1/(γ - cos(2π*exp(q)))

function generate_logfrq_pitch_comb(γ::Real, K::Int; frq_per_octave::Number = 1000)
    frqs = SpeechBox.logfreq_array(fmin = 0.5, fmax = K + 0.5, frq_per_octave = frq_per_octave)
    q = log.(frqs)
    h = zeros(size(q))
    for i = 1:length(q)
        h[i] = pitch_estimation_filter(q[i],γ,K)
    end
    β = sum(h)/length(h)
    h = h .- β
    zindx = argmin(abs.(q))
    return h, zindx
end

function generate_pitch_estimate(x::Array{Float64,1},fs::Number)
    frq_per_octave = 1000
    γ = 1.8
    K = 5
    frqs = SpeechBox.logfreq_array(fmin = 1, fmax = fs/2, frq_per_octave = frq_per_octave)
    pd = periodogram(x,fs,frqs)
    h,z = generate_logfrq_pitch_comb(γ,K, frq_per_octave = frq_per_octave)
    padded_pd = zeros(length(pd.components)+length(h)-1)
    padded_pd[z:z+length(pd.components)-1] = log.(pd.components)
    y = zeros(size(pd.components))
    for i=1:length(comb_output)
        y[i] = dot(padded_pd[i:i+length(h)-1],h)
    end
    return y
end

function LTASS(frqs::Vector{T}) where T<:Real
    f = [63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000]
    m = [38.6, 43.5, 54.4, 57.7, 56.8, 60.2, 60.3, 59.0, 62.1, 62.1, 60.5, 56.8, 53.7, 53.0, 52.0, 48.7, 48.1, 46.8, 45.6, 44.5, 44.3, 43.7, 43.4, 41.3, 40.7]
    itp = interpolate((f,),m,Gridded(Linear()))
    ltass_mag = zeros(size(frqs))
    for i=1:length(frqs)
        ltass_mag[i] = itp(frqs[i])
    end
    return ltass_mag
end

function LTASS(fmax::Number, frq_per_octave::Number)
    f = [63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000]
    m = [38.6, 43.5, 54.4, 57.7, 56.8, 60.2, 60.3, 59.0, 62.1, 62.1, 60.5, 56.8, 53.7, 53.0, 52.0, 48.7, 48.1, 46.8, 45.6, 44.5, 44.3, 43.7, 43.4, 41.3, 40.7]
    itp = interpolate((f,),m,Gridded(Linear()))
    frqs = logfreq_array(fmin = f[1], fmax = min(fmax,f[end]), frq_per_octave = 1000)
    ltass_mag = zeros(size(frqs))
    for i=1:length(frqs)
        ltass_mag[i] = itp(frqs[i])
    end
    return ltass_mag, frqs
end



function moving_average(pd::timefreq; time_window::Int, freq_window::Int)
    tfilt = DSP.Filters.PolynomialRatio(1/time_window*ones(time_window),[1])
    ffilt = DSP.Filters.PolynomialRatio(1/freq_window*ones(freq_window),[1])
    smoothed_components = zeros(size(pd.components))
    nfrqs,nframes = size(pd.components)
    time_pad_len = Int((time_window-1)/2)
    freq_pad_len = Int((freq_window-1)/2)
    for i=1:nfrqs
        ip = symmetric_pad(pd.components[i,:],time_pad_len)
        op = DSP.filt(tfilt,ip)
        smoothed_components[i,:] = unpad_vector(op,time_pad_len)
    end
    for i=1:nframes
        ip = symmetric_pad(smoothed_components[:,i],freq_pad_len)
        op = DSP.filt(ffilt,ip)
        smoothed_components[:,i] = unpad_vector(op,freq_pad_len)
    end
    return timefreq(pd.frames,smoothed_components,pd.frqs)
end


function LTASS_normalise(pd::timefreq)
    time_window = 11
    freq_window = 11
    frq_per_octave = 1000
    ltass_mag = LTASS(pd.frqs)
    smpd = moving_average(pd,time_window = time_window, freq_window = freq_window)
    norm_components = zeros(size(pd.components))
    for i=1:size(pd.components,2)
        norm_components[:,i] = pd.components[:,i].*ltass_mag./smpd.components[:,i]
    end
    return timefreq(pd.frames,norm_components,pd.frqs)
end
