#Speech waveform (time-domain) object
struct speech_waveform
    x::Array{Float,1}
    fs::Float
end

#Constructor for speech_waveform object, taking a Vector or Matrix and sampling frequency as input
function speech_waveform(x::Array,fs::Number)
    fs = convert(Float,fs)
    if(eltype(x)!=Float)
        error("Input signal is not Float, check if conversion is required")
    end
    if(ndims(x)>2)
        error("Signal has more than 2 dimensions, interpretation is not known")
    end
    if(ndims(x)==2)
        println("Input is a Matrix, assuming more than one channel and taking the first one")
        if(size(x,1)>size(x,2))
            signal = x[:,1]
        else
            signal = x[1,:]
        end
    else
        signal = x
    end
    return speech_waveform(signal,fs)
end


#Struct/Object defining framed signal using signal defined as a 1D array
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
function framed_signal(x::Array{<:AbstractFloat},fs::Number,win_dur::AbstractFloat=0.02,win_shift::AbstractFloat=0.01)
    signal = speech_waveform(x,fs)
    return framed_signal(signal,win_dur,win_shift)
end

#Union types for real or complex vector/matrix to denote spectral and spectro-temporal representations
RorC_Vector = Union{Array{Float,1},Array{Complex{Float},1}}
RorC_Matrix = Union{Array{Float,2},Array{Complex{Float},2}}

#Spectrum object to hold any form of frequency components of a signal
struct spectrum
    signal::speech_waveform
    components::RorC_Vector
    frqs::Array{Float,1}
    title::AbstractString
end

spectrum(signal::speech_waveform,components::RorC_Vector,frqs::Vector{Float}) = spectrum(signal,components,frqs,"")
