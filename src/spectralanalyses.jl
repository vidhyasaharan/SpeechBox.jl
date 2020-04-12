#Function to estimate Fourier spectrum of signal using DFT
function dftspec(signal::speech_waveform;wtype::String="hanning")
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    x = signal.x
    flen = length(x)
    win = window(flen;wtype=wtype)
    cmplx_spectrum = rfft(x.*win)
    nfrqs = length(cmplx_spectrum)
    frqs = range(0, signal.fs, length = nfrqs)
    return spectrum(signal,cmplx_spectrum,frqs,"Complex Fourier Spectrum")
end

function dftspec(x::Array{Float,1},fs::Number=1.0;wtype::String="hanning")
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    signal = speech_waveform(x,fs)
    return dftspec(signal;wtype=wtype)
end

function magspec(signal::speech_waveform;wtype::String="hanning")
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    cspec = dftspec(signal;wtype=wtype)
    return spectrum(signal,abs.(cspec.components),cspec.frqs,"DFT Magnitude Spectrum")
end

function magspec(x::Array{Float,1},fs::Float=1.0;wtype::String="hanning")
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    signal = speech_waveform(x,fs)
    return magspec(signal;wtype=wtype)
end

#Spectrogram estimated from framed_signal object input (core method for later verions)
function specgram(sig_frames::framed_signal;wtype::String="hanning")
    flen = sig_frames.frame_length
    nfft = nextfastfft(flen) #Get optimal number of points (larger than frame length) for FFT
    nframes = sig_frames.num_signal_frames

    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    win = window(flen;wtype=wtype)

    buf = zeros(nfft); #Buffer for operating on one frame (length is equal or larger than frame length)
    rfp = plan_rfft(buf); #Real valued FFT operator (gives only positive frequencies)
    nrfft = length(rfp*buf); #Number of FFT coefficeints
    mspec = zeros(nrfft,nframes); #Buffer for spectrogram values
    for i=1:nframes
        frame = extract_frame(sig_frames,i)
        buf[1:1:flen] = win.*frame; #Apply window and THEN store in buffer
        mspec[:,i] = abs.(rfp*buf); #Magnitude spectrum
    end
    return mspec + (eps()*ones(size(mspec))) #Add a tiny floor to spectrogram to avoid potential zero values - in case log spectrogram is required later.
end

#Spectrogram wrapper for Array{AbstactFloat} input
function specgram(x::Array{<:AbstractFloat},fs::Number;win_dur::Float=0.02,win_shift::Float=0.01,wtype::String="hanning")
    sig_frames = framed_signal(x,fs,win_dur,win_shift) #Obtain signal frames object
    return specgram(sig_frames;wtype=wtype)
end




#Periodogram estimated at provided frequncies - estimated by projecting onto complex exponentials and taking the square of the absolute value
function periodogram(x::Array{Float,1},fs::Number,frqs::Array{T,1};wtype::String="hanning") where T<:Number
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    flen = length(x)
    win = window(flen;wtype=wtype)
    proj_matrix = cexp_proj_matrix(frqs,fs,flen)
    proj = proj_matrix*(x.*win)
    return abs2.(proj)
end


#Wrapper function for periodogram over logarithmically spaced frequencies
function periodogram(x::Array{Float,1},fs::Number;wtype::String="hanning",fmin::Number=10,fmax::Number=fs/2)
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram(x,fs,frqs;wtype=wtype)
end
