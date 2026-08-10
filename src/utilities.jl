#Speech-specific signal conditioning. Generic signal utilities (windows, resampling,
#frequency grids, signal generators) live in TimeFrequencyAnalysis.jl.

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
