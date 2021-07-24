# Function that pulls out one frame as an array from framed_signal object
"""
    extract_frame(x::framed_signal, i)

Extract frame number `i` from the the framed\\_signal object `x`

"""
extract_frame(x::framed_signal,i::Int) = collect(view_frame(x, i))
    


# Function that creates a view of one frame as a subarray from framed_signal object
"""
    view_frame(x::framed_signal, i)

Create of view of frame number `i` from the the framed\\_signal object `x` (Frames involving zero padding of the signal will be arrays not views).

"""
function view_frame(x::framed_signal,i::Int)
    frame_length = x.frame_length
    frame_shift = x.frame_shift
    sindx = (i-1)*frame_shift + 1 #Identify start index of desired frame
    eindx = sindx + frame_length - 1 #Identify end index of desired frame
    if(i<=x.num_signal_frames) #Checking to see end of frame is within bounds of defined signal
        frame = @view x.signal.x[sindx:eindx]
    else #Zero padding to return full frame if signal ends midway through the frame
        frame = zeros(Float,frame_length)
        lindx = length(x.signal.x)
        sig_len = lindx - sindx + 1
        frame[1:sig_len] = @view x.signal.x[sindx:lindx]
    end
    return frame
end



#Compute the number of complete frames in signal (without zero padding) given frame size and frame shift
"""
    number_signal_frames(s::speech_waveform, frame_size, frame_shift)
    number_signal_frames(s::speech_waveform, frame_dur, frame_shift_dur)

Compute the number of frames without zero padding or extension of speech waveform `s` given `frame_size` and `frame_shift` in number of samples or `frame_dur` and `frame_shift_dur` in seconds
"""
number_signal_frames(s::speech_waveform, frame_size::Int, frame_shift::Int) = 1+ Int(floor((length(s.x)-frame_size)/frame_shift))

function number_signal_frames(s::speech_waveform, frame_dur::Float, frame_shift_dur::Float)
    frame_size = time2nsamples(frame_dur,s.fs)
    frame_shift = time2nsamples(frame_shift_dur,s.fs)
    return number_signal_frames(s, frame_size, frame_shift)
end


"""
    frame_energy(sig_frames::framed_signal[; normalised = false])

Estimate energy in each frame of `sig_frames` as sum of squares of samples (optionally return normalised energy such that max is 1 if `normalised = true`)

"""
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
