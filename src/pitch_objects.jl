#Object pairing a pitch contour with a time-frequency representation, so that the plot
#recipe in plot_recipes.jl can overlay the contour on a spectrogram.

"""
    pitch_timefreq(pitch, msp::timefreq)
    pitch_timefreq(pitch, frames::framed_signal)

Pair a `pitch` contour (one f₀ value in Hz per frame, `NaN` for unvoiced frames) with the
time-frequency representation `msp` for plotting. `pindx` holds the row index of `msp`
closest to each pitch value (`NaN` where unvoiced). Given a [`framed_signal`](@ref) instead
of a `timefreq`, the dB magnitude spectrogram of the frames is computed and used.

Returns `nothing` if the length of `pitch` does not match the number of frames.
"""
struct pitch_timefreq{T<:AbstractFloat, S<:Union{T,Complex{T}}}
    pitch::Vector{T}
    pindx::Vector{T}
    msp::timefreq{T,S}
end

function pitch_timefreq(pitch::AbstractVector{<:Real}, msp::timefreq{T}) where {T<:AbstractFloat}
    if(length(pitch)==size(msp.components,2))
        pindex = Vector{T}(undef,length(pitch))
        fill!(pindex,T(NaN))
        for i ∈ eachindex(pitch)
            if(~isnan(pitch[i]))
                pindex[i] = T(frqindex(pitch[i],msp.frqs))
            end
        end
        return pitch_timefreq(convert(Vector{T},pitch), pindex, msp)
    end
end

function pitch_timefreq(pitch::AbstractVector{<:Real}, frames::framed_signal)
    if(length(pitch)==frames.num_signal_frames)
        msp = specgram(frames)
        return pitch_timefreq(pitch, amp2db(msp))
    end
end
