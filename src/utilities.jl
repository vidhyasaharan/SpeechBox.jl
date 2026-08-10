#Pre-emphasis with 1 - 0.95z^-1
"""
    preemphasis(x::AbstractVector)
    preemphasis(s::speech_waveform)

Pre-emphasise speech signal stored in array `x` or given by `speech_waveform` object `s` using high pass filter: ``H(z) = 1 - 0.95z^-1``
"""
function preemphasis(x::AbstractVector{T}) where {T<:AbstractFloat}
    y = Vector{T}(undef,length(x))
    α = T(0.95)
    y[1] = T(0.05)*x[1]
    @inbounds @simd for i ∈ 2:length(x)
        y[i] = x[i] - α*x[i-1]
    end
    return y
end

preemphasis(s::speech_waveform) = speech_waveform(preemphasis(s.x),s.fs)



#Generate window function of given length
"""
    window(len [; wtype = "hanning"])
    window(T, len [; wtype = "hanning"])

Generate a window vector of length `len` and type `wtype`, with element type `T<:AbstractFloat` (default is `Float64`). Default window type is the Hann window. Options for wtype are:

### Implmented window types (`wtype`)
- "hanning" : Hann window [Default]
- "hamming" : Hamming window
- "rect" : Rectangular window
"""
function window(::Type{T}, flen::Int; wtype::String="hanning") where {T<:AbstractFloat}
    if(wtype=="rect")
        win = ones(T,flen)
    elseif(wtype=="hamming")
        win = convert(Vector{T},hamming(flen))
    elseif(wtype=="hanning")
        win = convert(Vector{T},hanning(flen))
    else
        println("Warning: window type not recognised - using Hann window")
        win = convert(Vector{T},hanning(flen))
    end
    return win
end

window(flen::Int;wtype::String="hanning") = window(Float64, flen; wtype)



#Resample signal in speech_waveform object (wrapper for resample from DSP.jl)
"""
    resample(s::speech_waveform, fs_new)

Resample the signal `s` to new sampling rate `fs_new` using the `resample` method from [`DSP.jl`](https://docs.juliadsp.org/stable/contents/)
"""
function resample(signal::speech_waveform, fs_new::Number)
    rx = DSP.Filters.resample(signal.x, fs_new/signal.fs)
    return speech_waveform(rx,fs_new)
end
