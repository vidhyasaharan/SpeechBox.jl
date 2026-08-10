struct spectral_comb end

#Frame based pitch estimation based on spectral comb
"""
    pitch(spectral_comb(), frames::framed_signal)

Compute the fundamental frequency (f₀) within each frame of a speech signal in the `frames` based on the position of peak of the spectral comb response. The `vad` function is used to determine voicing in each frame and `NaN` is returned for unvoiced frames.
"""
function pitch(::spectral_comb, frames::framed_signal{T}) where {T<:AbstractFloat}
    max_f₀ = 600.0  #Assume highest f₀ is 200Hz
    comb_resp, frqs = xcorr_spectral_comb(frames)
    mindx = frqindex(max_f₀,frqs)
    p = Vector{T}(undef,frames.num_signal_frames)
    fill!(p,T(NaN))
    v = vad(energy_threshold(), frames, threshold=0.02)
    for i ∈ axes(comb_resp,2)
        pindx = argmax(comb_resp[1:mindx,i])
        if(v[i]==1)
            p[i] = frqs[pindx]
        end
    end
    return p
end