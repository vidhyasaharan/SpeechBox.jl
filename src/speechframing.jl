#Struct/Object defining framed signal using signal defined as SampleBuf from SampledSignals toolbox
struct framed_signal
    x::SampleBuf
    frame_length::Int
    frame_overlap::Int
    num_signal_frames::Int
    num_frames::Int
end

#Constructor function - sets up framed_signal object given an input signal and window parameters
function framed_signal(x::SampleBuf,win_dur::Float=0.02,win_overlap::Float=0.01)
    fs = samplerate(x)
    frame_length = Int(round(win_dur*fs))
    frame_overlap = Int(round(win_overlap*fs))
    num_signal_frames = Int(floor(1 + (nframes(x)-frame_length)/frame_overlap))
    num_frames = Int(ceil(nframes(x)/frame_overlap))
    return framed_signal(x,frame_length,frame_overlap,num_signal_frames,num_frames)
end

# Function that pulls out one frame as an array from framed_signal object
function extract_frame(x::framed_signal,i::Int)
    frame_length = x.frame_length
    frame_overlap = x.frame_overlap
    sindx = (i-1)*frame_overlap + 1 #Identify start index of desired frame
    eindx = sindx + frame_length - 1 #Identify end index of desired frame
    if(i<=x.num_signal_frames) #Checking to see end of frame is within bounds of defined signal
        frame = float(collect(x.x[sindx:eindx]))
    else #Zero padding to return full frame if signal ends midway through the frame
        frame = zeros(frame_length)
        lindx = nframes(x.x)
        sig_len = lindx - sindx + 1
        frame[1:sig_len] = float(collect(x.x[sindx:lindx]))
    end
    return frame
end

#Estimate energy in each frame as sum of squares of samples (optionally return normalised energy such that max is 1)
function frame_energy(sig_frames::framed_signal; normalised = false)
    numframes = sig_frames.num_frames
    energy = zeros(Float,numframes) #Initialise array of estimated energy values (one per frames)
    for i=1:numframes
        frame = extract_frame(sig_frames,i) #extract frame from framed signal object
        energy[i] = sum(abs2,frame) #Estimate energy as sum of squares
    end
    if normalised #if normalised flag is true, normalise energy such that maximum is 1
        return energy/maximum(energy)
    else
        return energy
    end
end
