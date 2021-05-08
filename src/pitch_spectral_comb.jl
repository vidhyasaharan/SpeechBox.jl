
#Frame based pitch estimation based on spectral comb
function pitch_spec_comb(frames::framed_signal)
    comb_resp, frqs = xcorr_spectral_comb(frames)
    p = Vector{Float}(undef,frames.num_signal_frames)
    fill!(p,NaN)
    v = vad(frames)
    for i ∈ axes(comb_resp,2)
        pindx = argmax(comb_resp[:,i])
        if(v[i]==1)
            p[i] = frqs[pindx]
        end
    end
    return p
end

