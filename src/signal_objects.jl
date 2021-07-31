#Speech waveform (time-domain) object
"""
    speech_waveform(x,fs)

Store the signal in array `x` with sampling rate `fs` as a speech waveform object. If `x` is a matrix, the longer dimension is chosen as the signal length and the first row or column of that lenght is stored.

"""
struct speech_waveform
    x::Vector{Float}
    fs::Float
end

#Constructor for speech_waveform object, taking a Vector or Matrix and sampling frequency as input
function speech_waveform(x::AbstractArray{<:AbstractFloat},fs::Number)
    fs = convert(Float,fs)::Float
    if(eltype(x)!=Float)
        error("Input signal is not Float, check if conversion is required")
    end
    if(ndims(x)>2)
        error("Signal has more than 2 dimensions, interpretation is not known")
    end
    if(ndims(x)==2)
        if((size(x,1)==1)||(size(x,2)==1))
            signal = x[:]
        else
            println("Input is a Matrix, assuming more than one channel and taking the first one")
            if(size(x,1)>size(x,2))
                signal = x[:,1]
            else
                signal = x[1,:]
            end
        end
    else
        signal = x
    end
    signal = convert(Vector{Float},signal)
    return speech_waveform(signal,fs)
end


#Struct/Object defining framed signal using signal defined as a 1D array
"""
    framed_signal(sig::speech_waveform, win_dur=0.02, win_shift=0.01)
    framed_signal(x, fs, win_dur=0.02, win_shift=0.01)

Store the signal in a speech\\_waveform object `sig` (or signal in array `x` with sampling rate `fs`) as a framed\\_signal object with frame duration `win_dur` (default = 0.02sec) and interval between start of consecutive frames `win_shift` (detault = 0.01sec)

"""
struct framed_signal
    signal::speech_waveform
    frame_length::Int
    frame_shift::Int
    num_signal_frames::Int
    num_frames::Int
end


#Consutrctor for framed_signal object given input speech_waveform object, window duration and window overlap
function framed_signal(signal::speech_waveform,win_dur::Number=0.02,win_shift::Number=0.01)
    fs = signal.fs
    x = signal.x
    frame_length = Int(round(win_dur*fs))
    frame_shift = Int(round(win_shift*fs))
    num_signal_frames = Int(floor(1 + (length(x)-frame_length)/frame_shift))
    num_frames = Int(ceil(length(x)/frame_shift))
    return framed_signal(signal,frame_length,frame_shift,num_signal_frames,num_frames)
end


#Constructor function - sets up framed_signal object given an input signal of one dimenionsal float array and window parameters
function framed_signal(x::Array{<:AbstractFloat},fs::Number,win_dur::Number=0.02,win_shift::Number=0.01)
    signal = speech_waveform(x,fs)
    return framed_signal(signal,win_dur,win_shift)
end

#Union types for real or complex vector/matrix to denote spectral and spectro-temporal representations
RorC_Vector = Union{Array{Float,1},Array{Complex{Float},1}}
RorC_Matrix = Union{Array{Float,2},Array{Complex{Float},2}}
RorC = Union{Float,Complex{Float}}

#Spectrum object to hold any form of frequency components of a signal
"""
    spectrum(signal, components, frqs, title)
    spectrum(signal, components, frqs)

Construct a spectrum object to hold the spectral elements in `components` corresponding to the speech\\_waveform `signal` with frequency indices `frqs` and stores `title`
"""
struct spectrum{T<:RorC}
    signal::speech_waveform
    components::Vector{T}
    frqs::Vector{Float}
    title::AbstractString
end

function spectrum(signal::speech_waveform, components::AbstractVector{<:RorC},frqs::AbstractVector{<:Real},title::AbstractString)
    frqs = convert(Vector{Float},frqs)
    return spectrum(signal,components,frqs,title)
end

spectrum(signal::speech_waveform,components::Vector{<:RorC},frqs::Vector{<:Real}) = spectrum(signal,components,frqs,"")

#Spectrum object to hold any form of spectro-temporal components of a signal
"""
    timefreq(sig, frames, components, frqs, time, title)
    timefreq(sig, components, frqs, time, title)
    timefreq(frames, components, frqs, time, title)
    timefreq(frames, components, frqs, title)
    timefreq(frames, components, frqs)

Construct a time-frequuency object to hold the time-frequency elements in `components` corresponding to the speech\\_waveform `signal` based on framed\\_signal `frames` with frequency indices `frqs`, time indices `time` and stores `title`
"""
struct timefreq
    signal::speech_waveform
    frames::Union{framed_signal, Nothing}
    components::RorC_Matrix
    frqs::Vector{Float}
    time::Union{Vector{Float}, Vector{Vector{Float}}} #Vector if time indices are consistent across components, Vector of Vector if components have different time
    title::Union{AbstractString, Nothing}
end

#Contructor functions for expected input combinations
timefreq(frames::framed_signal,
        components::RorC_Matrix,
        frqs::Vector{Float},
        time::Union{Vector{Float}, Vector{Vector{Float}}},
        title::Union{AbstractString, Nothing}) = timefreq(frames.signal,frames,components,frqs,time,title)

timefreq(signal::speech_waveform,
        components::RorC_Matrix,
        frqs::Vector{Float},
        time::Union{Vector{Float}, Vector{Vector{Float}}},
        title::Union{AbstractString, Nothing}) = timefreq(signal,nothing,components,frqs,time,title)

function timefreq(frames::framed_signal, components::RorC_Matrix, frqs::Vector{Float}, title::Union{AbstractString, Nothing})
    frame_shift = frames.frame_shift
    fs = frames.signal.fs
    frame_length = frames.frame_length
    num_frames = size(components,2)

    time_shift = frame_shift/fs
    start_time = frame_length/(2*fs)
    t = collect(range(start_time, step = time_shift, length = num_frames))
    return timefreq(frames.signal,frames,components,frqs,t,title)
end

timefreq(frames::framed_signal, components::RorC_Matrix, frqs::Vector{Float}) =timefreq(frames,components,frqs,nothing)



#Object to hold pitch contour and spectrogram to allow plotting

struct pitch_timefreq
    pitch::Vector{Float}
    pindx::Vector{Float}
    msp::timefreq
end

function pitch_timefreq(pitch::Vector{Float}, msp::timefreq)
    if(length(pitch)==size(msp.components,2))
        pindex = Vector{Float}(undef,length(pitch))
        fill!(pindex,NaN)
        for i ∈ eachindex(pitch)
            if(~isnan(pitch[i]))
                pindex[i] = Float(frqindex(pitch[i],msp.frqs))
            end
        end
        return pitch_timefreq(pitch, pindex, msp)
    end
end

function pitch_timefreq(pitch::Vector{Float}, frames::framed_signal)
    if(length(pitch)==frames.num_signal_frames)
        msp = specgram(frames)
        return pitch_timefreq(pitch, amp2db(msp))
    end
end