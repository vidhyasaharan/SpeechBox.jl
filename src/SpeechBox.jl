module SpeechBox


using Reexport
using DSP
using FFTW
using StaticArrays

@reexport using SampledSignals

export framed_signal, extract_frame, specgram, melbankm, melfcc, Float, vad_energy_threshold

const Float = Float64 #Set to Float32 to fore 32bit floating point operations - NOT YET IMPLEMENTED or TESTED

struct framed_signal
    x::SampleBuf
    frame_length::Int
    frame_overlap::Int
    num_signal_frames::Int
    num_frames::Int
end

function framed_signal(x::SampleBuf,win_dur=0.02,win_overlap=0.01)
    fs = samplerate(x)
    frame_length = Int(round(win_dur*fs))
    frame_overlap = Int(round(win_overlap*fs))
    num_signal_frames = Int(floor(1 + (nframes(x)-frame_length)/frame_overlap))
    num_frames = Int(ceil(nframes(x)/frame_overlap))
    return framed_signal(x,frame_length,frame_overlap,num_signal_frames,num_frames)
end

function extract_frame(x::framed_signal,i::Int)
    frame_length = x.frame_length
    frame_overlap = x.frame_overlap
    sindx = (i-1)*frame_overlap + 1
    eindx = sindx + frame_length - 1
    if(i<=x.num_signal_frames)
        frame = float(collect(x.x[sindx:eindx]))
    else
        frame = zeros(frame_length)
        lindx = nframes(x.x)
        sig_len = lindx - sindx + 1
        frame[1:sig_len] = float(collect(x.x[sindx:lindx]))
    end
    return frame
end

function vad_energy_threshold(sig_frames::framed_signal,energy_thr::Float)
    nframes = sig_frames.num_frames
    energy = zeros(Float,nframes)
    for i=1:nframes
        frame = extract_frame(sig_frames,i)
        energy[i] = rms(frame);
    end
    max_energy = maximum(energy)
    abs_thr = max_energy*energy_thr
    v_indx = zeros(Int,nframes)
    for i in findall(energy.>abs_thr)
        setindex!(v_indx,1,i)
    end
    return v_indx
end

function specgram(x::SampleBuf;win_dur=0.02,win_overlap=0.01,wtype="hanning")
    sig_frames = framed_signal(x,win_dur,win_overlap) #Obtain signal frames object
    flen = sig_frames.frame_length
    nfft = nextfastfft(flen) #Get optimal number of points (larger than frame length) for FFT
    nframes = sig_frames.num_signal_frames

    #Choose window - options are rectangle, hamming or hanning (function default is hanning)
    if(wtype=="rect")
        win = ones(flen);
    elseif(wtype=="hamming")
        win = hamming(flen);
    else
        win = hanning(flen)
    end

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

function melfcc(x::SampleBuf;ncoef=13,nfilt=17,win_dur=0.02,win_overlap=0.01)
    frames = framed_signal(x,win_dur,win_overlap)
    fs = samplerate(x)
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


function melbankm(fs,npts;nfilt=17)
    mspace = SVector{nfilt}(range(0,stop=frq2mel(fs/2),length = nfilt)) #equally spaced points in mel scale
    fspace = SVector{nfilt}(mel2frq.(mspace)) #equal mel spaced points mapped back to Hz
    cindx = SVector{nfilt}(Int.(round.(fspace*(npts-1)/(fs/2))+1)) #Filter centre indices (first at 0, final at Fs/2)
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

frq2mel(frq) = log(1+frq/700)*1127.01048;

mel2frq(mel) = 700*(exp(mel/1127.01048)-1);


end # module
