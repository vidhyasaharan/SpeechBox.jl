
#Autocorrelation function at time index `i` and lag `k` (in terms of samples)
"""
    acf(sig, i, k; win_size)
    acf(sig, t, lags; win_dur)

Computes the autocorrelation function of input signal `sig` provided as a `speech_waveform` object at sample `i` (or at time `t` in secs) at one of more lag index/indices `k` (or at time lags `lags` in secs) given a window size of `win_size` samples (or a window duration `win_dur` given in seconds)
"""
acf(s::speech_waveform, i::Int, k::Int; win_size::Int) = dotavx(s.x[i:i+win_size-k-1], s.x[i+k:i+win_size-1])

function acf(s::speech_waveform, t::Real, lag::Real; win_dur::Real)
    i = time2nsamples(t,s.fs)
    k = time2nsamples(lag,s.fs)
    win_size = time2nsamples(win_dur,s.fs)
    return acf(s,i,k;win_size)
end

acf(s::speech_waveform, i::Int, k::AbstractVector{Int}; win_size::Int) = map(x->acf(s,i,x;win_size),k)

acf(s::speech_waveform, t::Real, lags::AbstractVector{<:Real}; win_dur::Real) = map(x->acf(s,t,x;win_dur),lags)



#Cross correlation function at time index `i` and lag `k`
"""
    nccf(sig, i, k; win_size)
    nccf(sig, t, lags; win_dur)

Computes the normalised cross correlation function (with a delayed version) of input signal `sig` provided as a `speech_waveform` object at sample `i` (or at time `t` in secs) at one of more lag index/indices `k` (or at time lags `lags` in secs) given a window size of `win_size` samples (or a window duration `win_dur` given in seconds)
"""
function nccf(s::speech_waveform, i::Int, k::Int; win_size::Int, nconst::Real = 0)
    s1 = s.x[i:i+win_size-1]
    s2 = s.x[i+k:i+k+win_size-1]
    mean_s = sum(s1)/length(s1)
    s1 .-= mean_s
    s2 .-= mean_s
    e1 = dotavx(s1)
    e2 = dotavx(s2)
    ccf = dotavx(s1, s2)
    return ccf/(sqrt(nconst + (e1*e2)))
end

function nccf(s::speech_waveform, t::Real, lag::Real; win_dur::Real, nconst::Float = 0.0)
    i = time2nsamples(t, s.fs)
    k = time2nsamples(lag, s.fs)
    win_size = time2nsamples(win_dur, s.fs)
    return nccf(s,i,k;win_size,nconst)
end

nccf(s::speech_waveform, i::Int, k::AbstractVector{Int}; win_size::Int, nconst::Float = 0.0) = map(x->nccf(s,i,x;win_size,nconst),k)

nccf(s::speech_waveform, t::Real, lags::AbstractVector{<:Real}; win_dur::Real, nconst::Float = 0.0) = map(x->nccf(s,t,x;win_dur,nconst),lags)


function nccf(s::speech_waveform; min_lag::Int, max_lag::Int, win_size::Int, win_shift::Int, nconst::Float = 0.0)
    lags = min_lag:max_lag
    nlags = length(lags)
    # nframes = 1 + Int(round((length(s.x) - (win_size+max_lag-1))/win_shift))
    nframes = number_signal_frames(s, win_size + max_lag, win_shift)
    cf = Matrix{Float}(undef,nlags,nframes)
    for j = 1:nframes
        for i ∈ eachindex(lags)
            sindx = (j-1)*win_shift + 1
            cf[i,j] = nccf(s, sindx, lags[i]; win_size, nconst)
        end
    end
    return cf
end