# Functions that fill a matrix of frames from a signal
"""
    enframe!(A, x)

Fills each column of matrix `A` with a frame of array `x`. Frame size is given by number of rows of `A` and frame shift is determined such that the final frame is as close as possible to end of `x`.
"""
function enframe!(y::AbstractMatrix{Float}, x::AbstractVector{Float})
    frame_size = size(y,1)
    num_frames = size(y,2)
    len = length(x)
    frame_shift = Int(floor((len - frame_size)/(num_frames-1)))
    for i ∈ 1:num_frames
        sindx = (i-1)*frame_shift + 1
        eindx = sindx + frame_size - 1
        y[:,i] = x[sindx:eindx]
    end
    return y
end


"""
    enframe(x, frame_size, frame_shift)

Generate a matrix with each column storing one frame of length `frame_size` from signal in array `x`, with consecutive frames separated by `frame_shift` samples.
"""
function enframe(x::AbstractVector{Float}, frame_size::Int, frame_shift::Int)
    num_frames = number_signal_frames(x, frame_size, frame_shift)
    y = Matrix{Float}(undef,frame_size,num_frames)
    enframe!(y, x)
    return y
end

"""
    enframe(x, fs, frame_dur, frame_shift_dur)

Generate a matrix with each column storing one frame of duration `frame_dur` in seconds from signal in array `x` with sampling rate `fs`, with consecutive frames separated by `frame_shift_dur` in seconds.
"""
function enframe(x::AbstractVector{Float}, fs::Real, frame_dur::Real, frame_shift_dur::Real)
    frame_size = time2nsamples(frame_dur,fs)
    frame_shift = time2nsamples(frame_shift_dur,fs)
    return enframe(x,frame_size, frame_shift)
end

"""
    enframe(signal::speech_waveform, frame_dur, frame_shift_dur)

Generate a matrix with each column storing one frame of duration `frame_dur` in seconds from speech\\_waveform `signal`, with consecutive frames separated by `frame_shift_dur` in seconds.
"""
enframe(s::speech_waveform, frame_dur::Real, frame_shift_dur::Real) = enframe(s.x, s.fs, frame_dur, frame_shift_dur)



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
    number_signal_frames(x, frame_size, frame_shift)
    number_signal_frames(x, fs, frame_dur, frame_shift_dur)
    number_signal_frames(s::speech_waveform, frame_size, frame_shift)
    number_signal_frames(s::speech_waveform, frame_dur, frame_shift_dur)

Compute the number of frames without zero padding or extension of array `x` or speech waveform `s` given `frame_size` and `frame_shift` in number of samples or `frame_dur` and `frame_shift_dur` in seconds
"""
number_signal_frames(x::AbstractVector{Float}, frame_size::Int, frame_shift::Int) = 1+ Int(floor((length(x)-frame_size)/frame_shift))

number_signal_frames(s::speech_waveform, frame_size::Int, frame_shift::Int) = number_signal_frames(s.x, frame_size, frame_shift)

function number_signal_frames(x::AbstractVector{Float}, fs::Real, frame_dur::Float, frame_shift_dur::Float)
    frame_size = time2nsamples(frame_dur, fs)
    frame_shift = time2nsamples(frame_shift_dur, fs)
    return number_signal_frames(x, frame_size, frame_shift)
end

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
