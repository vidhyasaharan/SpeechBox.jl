#Function to estimate Fourier spectrum of signal using DFT
"""
    dftspec(signal::speech_waveform[; wtype="hanning"])
    dftspec(x, fs=1.0[; wtype="hanning"])

Compute the complex DFT spectrum of speech\\_waveform `signal` (or signal in array `x` with sampling rate `fs`) using a window of type `wtype` (default = hanning window)

"""
function dftspec(signal::speech_waveform;wtype::String="hanning")
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    cmplx_spectrum,frqs = dftspec_components(signal.x, signal.fs; wtype)
    return spectrum(signal,cmplx_spectrum,frqs,"Complex Fourier Spectrum")
end

function dftspec(x::AbstractVector{Float},fs::Real=1.0;wtype::String="hanning")
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    signal = speech_waveform(x,fs)
    return dftspec(signal;wtype)
end


function dftspec_components(x::AbstractVector{Float}, fs::Real = 1.0; wtype::String="hanning")
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    flen = length(x)
    win = window(flen;wtype=wtype)
    cmplx_spectrum = rfft(x.*win)
    nfrqs = length(cmplx_spectrum)
    frqs = linfreq_array(fmin = 0, fmax = fs/2; nfrqs)
    # frqs = collect(range(0, fs/2, length = nfrqs))
    return cmplx_spectrum, frqs
end


"""
    magspec(signal::speech_waveform[; wtype="hanning"])
    magspec(x, fs=1.0[; wtype="hanning"])

Compute the DFT magnitude spectrum of speech\\_waveform `signal` (or signal in array `x` with sampling rate `fs`) using a window of type `wtype` (default = hanning window)

"""
function magspec(signal::speech_waveform;wtype::String="hanning")
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    cspec,frqs = dftspec_components(signal.x, signal.fs; wtype)
    return spectrum(signal,abs.(cspec),frqs,"DFT Magnitude Spectrum")
end

function magspec(x::AbstractVector{Float},fs::Real=1.0;wtype::String="hanning")
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    signal = speech_waveform(x,fs)
    return magspec(signal;wtype=wtype)
end

#Spectrogram estimated from framed_signal object input (core method for later verions)
"""
    specgram(frames[; wtype="hanning"])
    specgram(x, fs[; win_dur=0.02[, win_shift=0.01[, wtype="hanning"]]])

Compute the DFT based spectrogram of signal in framed\\_signal `frames` (or signal in array `x` with sampling rate `fs` using frames of duration `win_dur` and interval between frames `win_shift`) using a window of type `wtype` (default = hanning window). 

"""
function specgram(sig_frames::framed_signal;wtype::String="hanning")
    mspec, frqs = specgram_components(sig_frames; wtype)
    return timefreq(sig_frames,mspec,frqs)
end

#Spectrogram wrapper for Array{AbstactFloat} input
function specgram(x::Array{<:AbstractFloat},fs::Number;win_dur::Float=0.02,win_shift::Float=0.01,wtype::String="hanning")
    sig_frames = framed_signal(x,fs,win_dur,win_shift) #Obtain signal frames object
    return specgram(sig_frames;wtype=wtype)
end



function specgram_components(sig_frames::framed_signal;wtype::String="hanning")
    flen = sig_frames.frame_length
    nfft = nextfastfft(flen) #Get optimal number of points (larger than frame length) for FFT
    nframes = sig_frames.num_signal_frames

    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    win = window(flen;wtype=wtype)

    buf = zeros(nfft); #Buffer for operating on one frame (length is equal or larger than frame length)
    rfp = plan_rfft(buf); #Real valued FFT operator (gives only positive frequencies)
    nrfft = length(rfp*buf); #Number of FFT coefficeints
    mspec = Matrix{Float}(undef,nrfft,nframes); #Buffer for spectrogram values
    for i=1:nframes
        frame = view_frame(sig_frames,i)
        buf[1:1:flen] = win.*frame; #Apply window and THEN store in buffer
        mspec[:,i] = abs.(rfp*buf); #Magnitude spectrum
    end
    dithered_mspec = mspec + (eps()*ones(size(mspec))) #Add a tiny floor to spectrogram to avoid potential zero values - in case log spectrogram is required later.
    frqs = linfreq_array(fmin = 0, fmax = sig_frames.signal.fs/2, nfrqs = nrfft)
    # frqs = convert.(Float,collect(range(0, sig_frames.signal.fs/2, length = nrfft)))
    return dithered_mspec, frqs
end



#Periodogram estimated at provided frequncies - estimated by projecting onto complex exponentials and taking the square of the absolute value
"""
    periodogram(x, fs, frqs[; wtype="hanning"])
    periodogram(x, fs [;wtype="hanning"[, fmin=10[, fmax=fs/2]]])
    periodogram(sig_frames::framed_signal, frqs[; wtype="hanning"])
    periodogram(sig_frames::framed_signal[; wtype="hanning"[, fmin=10[, fmax=sig_frames.signal.fs/2]]])

Compute the periodogram of a signal in array `x` with sampling frequency `fs` at frequencies specified in `frqs` or frequencies equally spaced on the log-scale between `fmin` and `fmax` as the L2 norm of the inner product between a complex exponential and `x`. When the input is a framed signal object `sig_frames`, the periodogram for each frame is computed.
"""
function periodogram(x::Array{Float,1},fs::Number,frqs::Array{T,1};wtype::String="hanning") where T<:Number
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    components = periodogram_components(x,fs,frqs;wtype)
    signal = speech_waveform(x,fs)
    return spectrum(signal,components,frqs,"Periodogram")
end

function periodogram(x::Array{Float,1},fs::Number;wtype::String="hanning",fmin::Number=10,fmax::Number=fs/2)
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram(x,fs,frqs;wtype=wtype)
end

function periodogram(sig_frames::framed_signal,frqs ;wtype::String="hanning")
    pspec = periodogram_components(sig_frames, frqs; wtype)
    return timefreq(sig_frames,pspec,frqs)
end

function periodogram(sig_frames::framed_signal; wtype::String="hanning", fmin=10,fmax=sig_frames.signal.fs/2)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram(sig_frames, frqs, wtype = wtype)
end



#Periodogram components estimated at provided frequncies - estimated by projecting onto complex exponentials and taking the square of the absolute value
function periodogram_components(x::Vector{Float},fs::Real,frqs::Vector{<:Real};wtype::String="hanning")
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    frqs = Float.(frqs)::Vector{Float}
    fs = Float(fs)::Float
    flen = length(x)
    win = window(flen;wtype=wtype)
    ip = x.*win
    nfrqs = length(frqs)
    proj = zeros(Float,nfrqs)
    for i ∈ eachindex(proj)
        proj[i] = abs2(dotavx(ip,cexp(frqs[i],fs,flen)))
    end
    return proj
end

function periodogram_components(x::Array{Float,1},fs::Real;wtype::String="hanning",fmin::Real=10,fmax::Real=fs/2)
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram_components(x,fs,frqs;wtype=wtype)
end

function periodogram_components(sig_frames::framed_signal, frqs ;wtype::String="hanning")
    frqs = Float.(frqs)
    flen = sig_frames.frame_length
    nframes = sig_frames.num_signal_frames
    fs = sig_frames.signal.fs
    pspec = Matrix{Float}(undef,length(frqs),nframes)
    nfrqs = length(frqs)
    win = window(flen;wtype=wtype)

    ce_array = collect(transpose(cexp_proj_matrix(frqs,fs,flen)))

    for i in 1:nframes
        ip = extract_frame(sig_frames,i)
        mulavx!(ip,win)
        @views for j in 1:nfrqs
            pspec[j,i] = abs2(dotavx(ip,ce_array[:,j]))
        end
    end
    return pspec
end

function periodogram_components(sig_frames::framed_signal; wtype::String="hanning", fmin=10,fmax=sig_frames.signal.fs/2)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram_components(sig_frames, frqs, wtype = wtype)
end