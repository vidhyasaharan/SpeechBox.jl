#Prob of voicing based bifference highest and next highest peak in comb response
function prob_voiced(frames::framed_signal)
    v1 = 10.0
    v2 = 2.0
    peak_diff = comb_resp_peak_diff(frames)
    peak_dist = dist_pitch_range(frames)
    p1 = exp.(-1.0./(v1*peak_diff))
    p2 = exp.(-v2.*peak_dist)
    p = p1.*p2
    return p
end



function comb_resp_peak_diff(frames::framed_signal)
    comb_resp, _ = xcorr_spectral_comb(frames)
    peak_diff = Vector{Float64}(undef,frames.num_signal_frames)
    for i ∈ eachindex(peak_diff)
        peaks = findpeaks_sorted(comb_resp[:,i]; num_peaks = 2)
        mag = peaks[2]
        peak_diff[i] = (mag[1] - mag[2])/mag[1]
    end
    return peak_diff
end


function dist_pitch_range(frames::framed_signal; min_f₀::Float64 = 40.0, max_f₀::Float64 = 500.0)
    comb_resp, frqs = xcorr_spectral_comb(frames)
    peak_dist = zeros(Float64,frames.num_signal_frames)
    for i ∈ eachindex(peak_dist)
        ind = argmax(comb_resp[:,i])
        pf = frqs[ind]
        ldist = min_f₀ - pf
        udist = pf - max_f₀
        dist = maximum([ldist, udist])
        peak_dist[i] = (dist + abs(dist))/2
    end
    return peak_dist.*(2/frames.signal.fs)
end


function comb_resp_dist(a::AbstractVector{Float64}, b::AbstractVector{Float64})
    anorm = 1/sum(abs2,a)
    bnorm = 1/sum(abs2,b)
    ā = anorm.*a
    b̄ = bnorm.*b
    return sum(abs2,ā-b̄)
end


function comb_resp_Δ(frames::framed_signal)
    comb_resp, _ = xcorr_spectral_comb(frames)
    Δ = Vector{Float64}(undef,frames.num_signal_frames-1)
    for i ∈ eachindex(Δ)
        Δ[i] = comb_resp_dist(comb_resp[:,i],comb_resp[:,i+1])
    end
    return Δ
end