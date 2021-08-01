#Pre-emphasis with 1 - 0.95z^-1
"""
    preemphasis(x::AbstractVector)
    preemphasis(s::speech_waveform)

Pre-emphasise speech signal stored in array `x` or given by `speech_waveform` object `s` using high pass filter: ``H(z) = 1 - 0.95z^-1``
"""
function preemphasis(x::AbstractVector{Float})
    y = Vector{Float}(undef,length(x))
    y[1] = 0.05*x[1]
    @turbo for i ∈ 2:length(x)
        y[i] = x[i] - 0.95*x[i-1]
    end
    return y
end

preemphasis(s::speech_waveform) = speech_waveform(preemphasis(s.x),s.fs)



#Generate window function of given length
"""
    window(len [; wtype = "hanning"])

Generate a window vector of length `len` and type `wtype`. Default window type is the Hann window. Options for wtype are:

### Implmented window types (`wtype`)
- "hanning" : Hann window [Default]
- "hamming" : Hamming window
- "rect" : Rectangular window
"""
function window(flen::Int;wtype::String="hanning")
    if(wtype=="rect")
        win = ones(flen)
    elseif(wtype=="hamming")
        win = hamming(flen)
    elseif(wtype=="hanning")
        win = hanning(flen)
    else
        println("Warning: window type not recognised - using Hann window")
        win = hanning(flen)
    end
    return win
end



#Resample signal in speech_waveform object (wrapper for resample from DSP.jl)
"""
    resample(s::speech_waveform, fs_new)

Resample the signal `s` to new sampling rate `fs_new` using the `resample` method from [`DSP.jl`](https://docs.juliadsp.org/stable/contents/)
"""
function resample(signal::speech_waveform, fs_new::Number)
    rx = DSP.Filters.resample(signal.x, fs_new/signal.fs)
    fs = convert(Float,fs_new)
    return speech_waveform(rx,fs)
end