
#Autocorrelation function at time index `i` and lag `k` (in terms of samples)
"""
    acf(sig, i, k; win_size)
    acf(sig, t, lags; win_dur)

Computes the autocorrelation function of input signal `sig` provided as a `speech_waveform` object at sample `i` (or at time `t` in secs) at one of more lag index/indices `k` (or at time lags `lags` in secs) given a window size of `win_size` samples (or a window duration `win_dur` given in seconds)
"""
acf(s::speech_waveform, i::Int, k::Int; win_size::Int) = dotavx(s.x[i:i+win_size-k-1], s.x[i+k:i+win_size-1])

function acf(s::speech_waveform, t::Real, lag::Real; win_dur::Real)
    i = timeindex(t,s.fs)
    k = timeindex(lag,s.fs)
    win_size = timeindex(win_dur,s.fs)
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
function nccf(s::speech_waveform, i::Int, k::Int; win_size::Int)
    s1 = s.x[i:i+win_size-1]
    s2 = s.x[i+k:i+k+win_size-1]
    e1 = dotavx(s1)
    e2 = dotavx(s2)
    ccf = dotavx(s1, s2)
    return ccf/(sqrt(e1*e2))
end

function nccf(s::speech_waveform, t::Real, lag::Real; win_dur::Real)
    i = timeindex(t, s.fs)
    k = timeindex(lag, s.fs)
    win_size = timeindex(win_dur, s.fs)
    return nccf(s,i,k;win_size)
end

nccf(s::speech_waveform, i::Int, k::AbstractVector{Int}; win_size::Int) = map(x->nccf(s,i,x;win_size),k)

nccf(s::speech_waveform, t::Real, lags::AbstractVector{<:Real}; win_dur::Real) = map(x->nccf(s,t,x;win_dur),lags)
