"""
    melfcc(frames; nfilt=17, ncoeff=13)
    melfcc(signal; nfilt=17, ncoeff=13, wind_dur=0.02, win_shift=0.01)
    melfcc(x, fs; nfilt=17, ncoeff=13, wind_dur=0.02, win_shift=0.01)

Compute a vector of mel frequency cepstral coefficients for each frame of the input signal
# Arguments
- `frames`: a `framed_signal` object that holds the input signal and details of frames
- `signal`: a `speech_waveform` object holding the input speech waveform
- `x`: an array holding input speech signal
- `fs`: sampling rate
- `nfilt`: number of triangular mel filters used in the MFCC extraction [default=17]
- `ncoeff`: number of desired MFCC coefficients (must be less than `nfilt`) including C₀ [default=13]
- `win_dur`: duration of speech frame in secs [default=0.02]
- `win_shift`: time interval between consecutive windows in secs [default=0.01]
"""
function melfcc(frames::framed_signal;ncoef=13,nfilt=17)
    # frames = framed_signal(x,win_dur,win_shift)
    fs = frames.signal.fs
    numframes = frames.num_frames;
    flen = frames.frame_length;
    nfft = nextfastfft(flen);
    win = hamming(flen);

    buf = zeros(nfft);
    fbuf = zeros(nfilt);
    rfp = plan_rfft(buf);
    dcp = plan_dct(fbuf);

    nrfft = length(rfp*buf);
    fbank = melbankm(fs,nrfft,nfilt=nfilt); #generate mel filterbank filters
    fbank = fbank.^2; #squaring triangular mel-filters to multiple with power spectrum

    mfcc = zeros(ncoef,numframes);

    for i=1:numframes
        buf[1:1:flen] = win.*extract_frame(frames,i);
        #fbuf = dcp*(fbank*(abs.(rfp*buf)));
        fbuf = fbank*(abs2.(rfp*buf) + eps()*ones(nrfft)); #multiplying with power spectrum (square of mag spectrum) and accumulating
        fbuf = dcp*(log.(fbuf));
        mfcc[:,i] = fbuf[1:ncoef];
    end
    return mfcc
end

function melfcc(signal::speech_waveform;ncoef=13,nfilt=17,win_dur=0.02,win_shift=0.01)
    frames = framed_signal(signal,win_dur,win_shift)
    return melfcc(frames,ncoef=ncoef,nfilt=nfilt)
end

function melfcc(x::Array{<:AbstractFloat},fs::AbstractFloat;ncoef=13,nfilt=17,win_dur=0.02,win_shift=0.01)
    frames = framed_signal(x,fs,win_dur,win_shift)
    return melfcc(frames,ncoef=ncoef,nfilt=nfilt)
end

"""
    melbankm(fs, npts[; nfilt = 17])

Compute a matrix of triangular mel filter responses to multiply with a DFT spectrum given a sampling frequency `fs`, number of DFT points `npts` and optionally the number of equally spaced mel filters `nfilt`
"""
function melbankm(fs,npts;nfilt=17)
    mspace = range(0,stop=frq2mel(fs/2),length = nfilt) #equally spaced points in mel scale
    fspace = mel2frq.(mspace) #equal mel spaced points mapped back to Hz
    cindx = Int.(round.(fspace*(npts-1)/(fs/2)).+1) #Filter centre indices (first at 0, final at Fs/2)
    fbank = zeros(nfilt,npts);
    #define triangular filters for 2 to nfilt-1
    for i=2:nfilt-1
        fbank[i,cindx[i-1]:cindx[i]] = range(0,stop=1,length=cindx[i]-cindx[i-1]+1) #upward slope of triangle
        fbank[i,cindx[i]:cindx[i+1]] = range(1,stop=0,length=cindx[i+1]-cindx[i]+1) #downward slope of triangle
    end
    #one sided triangular filters for first and last filters
    fbank[1,cindx[1]:cindx[2]] = range(1,stop=0,length=cindx[2]-cindx[1]+1)
    fbank[nfilt,cindx[nfilt-1]:cindx[nfilt]] = range(0,stop=1,length=cindx[nfilt]-cindx[nfilt-1]+1)
    return fbank
end


"""
    frq2mel(f)

Convert frequency `f` from Hz to mel scale
"""
frq2mel(frq) = log(1+frq/700)*1127.01048;


"""
    mel2frq(m)

Convert frequency `m` from mel scale to Hz
"""
mel2frq(mel) = 700*(exp(mel/1127.01048)-1);
