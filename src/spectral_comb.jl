#Function implements point value of pitch comb filter
pitch_estimation_filter(q::T, γ::T) where {T} =  1/(γ - cos(2π*exp(q)))

#Function to generate comb filter for pitch estimation in log frequency domain
function generate_logfrq_pitch_comb(;γ::Real = 1.8, K::Int = 5,  frq_per_octave::Number = 1000)
    frqs = logfreq_array(fmin = 0.5, fmax = K + 0.5, frq_per_octave = frq_per_octave)
    q = log.(frqs)
    h = zeros(size(q))
    for i = 1:length(q)
        h[i] = pitch_estimation_filter(q[i],γ)
    end
    β = sum(h)/length(h)
    h = h .- β
    zindx = argmin(abs.(q))
    return h, zindx
end

#Function to generate pitch estimates for one frame by convolving with comb filter in log freq domain
function generate_pitch_estimate(x::Vector{Float},fs::Number)
    frq_per_octave = 1000
    # γ = 1.8
    # K = 5
    frqs = logfreq_array(fmin = 1, fmax = fs/2, frq_per_octave = frq_per_octave)
    pd = periodogram(x,fs,frqs)
    h,z = generate_logfrq_pitch_comb(;frq_per_octave)
    y = SpeechBox.xcorr(log.(pd.components),h,z)
    # padded_pd = zeros(length(pd.components)+length(h)-1)
    # padded_pd[z:z+length(pd.components)-1] = log.(pd.components)
    # y = zeros(size(pd.components))
    # for i ∈ eachindex(y)
    #     y[i] = SpeechBox.dotavx(padded_pd[i:i+length(h)-1],h)
    # end
    f₀ = frqs[argmax(y)]
    return f₀
end


#Function to generate a sequence of pitch estimates, one per frame, using a comb filter in the log freq domain
function generate_logfrq_pitch_comb(frames::framed_signal)
    frq_per_octave = 1000
    frqs = logfreq_array(fmin = 1, fmax = fs/2, frq_per_octave = frq_per_octave)
    pd = periodogram(frames,frqs)
end



