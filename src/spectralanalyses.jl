function magspec(x::Array{Float,1},fs::Float=1.0;wtype::String="hanning")
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    flen = length(x)
    win = window(flen;wtype=wtype)
    return abs.(rfft(x.*win))
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

#Spectrogram wrapper for SampleBuf input
# function specgram(x::SampleBuf;win_dur::Float=0.02,win_overlap::Float=0.01,wtype::String="hanning")
#     sig_frames = framed_signal(x,win_dur,win_overlap) #Obtain signal frames object
#     return specgram(sig_frames;wtype=wtype)
# end

#Spectrogram wrapper for Array{AbstactFloat} input
function specgram(x::Array{<:AbstractFloat},fs::AbstractFloat;win_dur::Float=0.02,win_overlap::Float=0.01,wtype::String="hanning")
    sig_frames = framed_signal(x,fs,win_dur,win_overlap) #Obtain signal frames object
    return specgram(sig_frames;wtype=wtype)
end


function periodogram(x::Array{Float,1},fs::Number;wtype::String="hanning",fmin::Number=10,fmax::Number=4000)
    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    flen = length(x)
    win = window(flen;wtype=wtype)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    periodogram_basis = periodogram_basis_matrix(frqs,fs,flen)
    proj = periodogram_basis*(x.*win)
    return abs2.(proj)
end

function periodogram_basis_matrix(frqs::Array{T,1},fs::Number,N::Int) where T<:Number
    nfrqs = length(frqs)
    periodogram_basis = zeros(Complex{Float},nfrqs,N)
    for i=1:nfrqs
        periodogram_basis[i,:] = cexp(frqs[i],fs,N)
    end
    return periodogram_basis
end
