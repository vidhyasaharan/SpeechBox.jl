#Estimate LPC order from sampling frequency
lpc_order(fs::Number) = Int(round(fs/1000))+2


#Compute LPC/AR coefficients for a discrete-time sequence
"""
    lpc(x, N)

Compute the Linear Predictive Coding (LPC) coefficients of order `N`, of a sequence `x`.
"""
function lpc(x::Array{Float,1},N::Int)
    a = ones(Float,N+1)
    temp,err = DSP.LPC.lpc(x,N)
    a[2:end] = temp
    # filter = DSP.Filters.PolynomialRatio([1],a)
    return a
end


#Compute the LPC/AR model magnitude response given a dicrete-time signal
function lpc_freqz(x::Array{Float,1}, fs::Number, N::Int = lpc_order(fs); frqs::Array{T,1} = linfreq_array(fmax = fs/2, nfrqs = length(x))) where T<:Number
    a = lpc(x,N)
    filter = DSP.Filters.PolynomialRatio([1],a)
    h = DSP.freqresp(filter, frqs * ((2pi) / fs))
    return h
end

#Compute the LPC/AR model magnitude response given a dicrete-time signal and store in spectrum object
"""
    lpc_response(x, fs[, N = round(fs/1000)+2 ]; frqs = linfreq_array(0, fs/2, length(x)))
    lpc_response(frames[, N = round(fs/1000)+2 ]; frqs = linfreq_array(0, fs/2, length(x)))

Compute the magnitude response of the Linear Predictive Coding (LPC) / Autoregressive (AR) filter model (of order `N`) of signal in array `x` with sampling rate `fs` at frequencies specified in `frqs`.
When the input is a framed signal `frames`, the magnitude response of the LPC/AR filter model in each frame is computed and concatenated to form an LPC spectrogram.
"""
function lpc_response(x::Array{Float,1}, fs::Number, N::Int = lpc_order(fs); frqs::Array{T,1} = linfreq_array(fmax = fs/2, nfrqs = length(x))) where T<:Number
    h = lpc_freqz(x, fs, N; frqs)
    return spectrum(speech_waveform(x,fs),abs.(h),frqs,"LPC/AR Model Magnitude Respose")
end


function lpc_response(sig_frames::framed_signal, N::Int = lpc_order(sig_frames.signal.fs); frqs::Array{T,1} = linfreq_array(fmax = sig_frames.signal.fs/2, nfrqs = sig_frames.frame_length)) where T<:Number
    nframes = sig_frames.num_signal_frames
    nfrqs = length(frqs)
    fs = sig_frames.signal.fs
    lpcspec = zeros(nfrqs,nframes)
    for i=1:nframes
        frame = extract_frame(sig_frames,i)
        lpcspec[:,i] = abs.(lpc_freqz(frame, fs, N; frqs))
    end
    return timefreq(sig_frames, lpcspec, frqs, "LPC/AR Spectrogram")
end



##LPC Utilities


#Generate allpole filter given pole frequencies, bandwidths and sampling frequency
function allpole(pf::Array{<:Number,1} = [1000, 1800, 2900, 3400, 5000, 6800], pbw::Array{<:Number,1} = [50, 120, 200, 300, 500, 800]; fs::Number = 16000)
    num_poles = length(pf)
    poles  = zeros(Complex{Float},2*num_poles)
    for i=1:num_poles
        r = exp(-pi*pbw[i]/fs)
        # println(r)
        poles[2i-1] = r*exp(1im*2*pi*pf[i]/fs)
        poles[2i] = conj(poles[2i-1])
    end
    return DSP.Filters.ZeroPoleGain([0],poles,1)
end

#Generate random resonance frequency and bandwidth for a vocal tract filter given a range for the resonance frequency
function rand_vocalfilter_resonance(fmin::Number = 0, fmax::Number = 8000)
    frequency = fmin + ((fmax-fmin)*rand())
    bw_min = 0.08*frequency
    bw_max = 0.2*frequency
    bandwidth = bw_min + ((bw_max-bw_min)*rand())
    return frequency, bandwidth
end

#Generate a 'random' allpole vocal tract filter model given a sampling frequency and number of resonances with the resonance frequencies and bandwidths chosen at random
function rand_allpole(fs::Number = 8000, num_res::Number = 10)
    fmax = fs/2
    fint = fmax/num_res
    frqs = zeros(Float,num_res)
    bws = zeros(Float,num_res)
    flo = 50
    fhi = 1.5*fint
    for i=1:num_res
        frqs[i],bws[i] = rand_vocalfilter_resonance(flo, fhi)
        fmid = i*fint + 0.5*fint
        flo = fmid-fint
        if(frqs[i]>flo)
            flo = frqs[i]
        end
        fhi = fmid+fint
    end
    ar_filt = allpole(frqs,bws; fs)
    return ar_filt, frqs, bws
end