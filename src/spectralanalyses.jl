"""
    dft(x [, <keyword argument>])

Compute the complex DFT coefficients of the array `x`.

### Keyword Arguments
- `ndft` : Number of DFT points [Default = `length(x)`]. If input `ndft` is less than the length of `x`, the number of DFT points is changed to `length(x)`
- `wtype` : Window type to use [Default = "hanning"]
"""
function dft(x::AbstractVector{Float}; ndft::Int = length(x), wtype::String="hanning")
    len = length(x)
    win = window(len;wtype)
    buflen = max(ndft,len)
    buf = zeros(Float, buflen)
    buf[1:len] = x.*win
    return rfft(buf)
end



"""
    magspec(signal [, <keywork arguments>])
    magspec(x, fs = 2π [, <keywork arguments>])

Compute the DFT magnitude spectrum of speech\\_waveform `signal` (or signal in array `x` with sampling rate `fs`). Returns `spectrum` object.

    magspec(comp(), x, fs [, <keyword arguments>])

Compute the DFT magnitude spectrum of signal in array `x`.

    magspec(comp(), x [, <keyword arguments>])

Compute the DFT magnitude spectrum of signal in array `x` with sampling rate `fs`. Returns array spectral components and array of frequencies

### Keyword Arguments
- `ndft` : Number of DFT points [Default value is length of signal array]
- `wtype` : Window type to use [Default is "hanning"]

"""
magspec(signal::speech_waveform;ndft::Int = length(signal.x), wtype::String="hanning") = magspec(signal.x, signal.fs; ndft, wtype)


function magspec(x::AbstractVector{Float},fs::Real=2π; ndft::Int = length(x), wtype::String="hanning")
    mspec,frqs = magspec(comp(), x, fs; ndft, wtype)
    return spectrum(speech_waveform(x,fs), mspec, frqs, "DFT Magnitude Spectrum")
end


function magspec(::comp, x::AbstractVector{Float},fs::Real; ndft::Int = length(x), wtype::String="hanning")
    mspec = magspec(comp(), x; ndft, wtype)
    frqs = rfftfreq(ndft, fs)
    return mspec,frqs
end

magspec(::comp, x::AbstractVector{Float}; ndft::Int = length(x), wtype::String = "hanning") = abs.(dft(x; ndft, wtype))


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

function specgram(::comp, x::AbstractVector{Float}, fs::Real; win_dur::Float=0.02, win_shift::Float=0.01, ndft::Int = nextfastfft(time2nsamples(win_dur,fs)), wtype::String="hanning")
    sig_frames = framed_signal(x,fs,win_dur,win_shift) #Obtain signal frames object
    return specgram(comp(), sig_frames;ndft, wtype)
end

function specgram(::comp, sig_frames::framed_signal; ndft::Int = nextfastfft(sig_frames.frame_length), wtype::String = "hanning")
    len = sig_frames.frame_length
    nfft = max(len,ndft) #Default value of ndft is the optimal number of points (larger than frame length) for FFT
    nframes = sig_frames.num_signal_frames #Spectrogram only covers frames witout zero padding original signal (i.e., may not include the last few samples of the signal)

    win = window(len;wtype=wtype) #Choose window - options are rectangle, hamming or hanning (function default is hanning)

    buf = zeros(nfft); #Buffer for operating on one frame (length is equal or larger than frame length)
    rfp = plan_rfft(buf); #Real valued FFT operator (gives only positive frequencies)
    nrfft = length(rfp*buf); #Number of FFT coefficeints
    mspec = Matrix{Float}(undef,nrfft,nframes); #Buffer for spectrogram values
    for i ∈ 1:nframes
        frame = view_frame(sig_frames,i)
        buf[1:1:len] = win.*frame; #Apply window and THEN store in buffer
        mspec[:,i] = abs.(rfp*buf); #Magnitude spectrum
    end

    #Add a tiny floor to spectrogram to avoid potential zero values - in case log spectrogram is required later.
    for i ∈ eachindex(mspec)
        mspec[i] += eps(Float)
    end
    return mspec
end

function specgram_components(sig_frames::framed_signal;wtype::String="hanning")
    len = sig_frames.frame_length
    nfft = nextfastfft(len) #Get optimal number of points (larger than frame length) for FFT
    nframes = sig_frames.num_signal_frames

    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    win = window(len;wtype=wtype)

    buf = zeros(nfft); #Buffer for operating on one frame (length is equal or larger than frame length)
    rfp = plan_rfft(buf); #Real valued FFT operator (gives only positive frequencies)
    nrfft = length(rfp*buf); #Number of FFT coefficeints
    mspec = Matrix{Float}(undef,nrfft,nframes); #Buffer for spectrogram values
    for i=1:nframes
        frame = view_frame(sig_frames,i)
        buf[1:1:len] = win.*frame; #Apply window and THEN store in buffer
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
    len = length(x)
    win = window(len;wtype=wtype)
    ip = x.*win
    nfrqs = length(frqs)
    proj = zeros(Float,nfrqs)
    for i ∈ eachindex(proj)
        proj[i] = abs2(dotavx(ip,cexp(frqs[i],fs,len)))
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
    len = sig_frames.frame_length
    nframes = sig_frames.num_signal_frames
    fs = sig_frames.signal.fs
    pspec = Matrix{Float}(undef,length(frqs),nframes)
    nfrqs = length(frqs)
    win = window(len;wtype=wtype)

    ce_array = collect(transpose(cexp_proj_matrix(frqs,fs,len)))

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