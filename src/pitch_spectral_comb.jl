
#Frame based pitch estimation based on spectral comb
function pitch_spec_comb(frames::framed_signal)
    comb_resp, frqs = xcorr_spectral_comb(frames)
    p = Vector{Float}(undef,frames.num_signal_frames)
    v = vad(frames)
    for i ∈ axes(comb_resp,2)
        pindx = argmax(comb_resp[:,i])
        p[i] = frqs[pindx]*v[i]
    end
    return p
end

