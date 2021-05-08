
#Frame based pitch estimation based on spectral comb
function pitch_spec_comb(frames::framed_signal)
    max_f₀ = 600.0  #Assume highest f₀ is 200Hz
    comb_resp, frqs = xcorr_spectral_comb(frames)
    mindx = frqindex(max_f₀,frqs)
    p = Vector{Float}(undef,frames.num_signal_frames)
    fill!(p,NaN)
    v = vad(frames)
    for i ∈ axes(comb_resp,2)
        pindx = argmax(comb_resp[1:mindx,i])
        if(v[i]==1)
            p[i] = frqs[pindx]
        end
    end
    return p
end

