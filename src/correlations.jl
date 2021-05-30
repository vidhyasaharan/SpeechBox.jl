
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
ccf(s::speech_waveform, i::Int, k::Int; win_size::Int) = dotavx(s.x[i:i+win_size-1], s.x[i+k:i+k+win_size-1])

function ccf(s::speech_waveform, t::Real, lag::Real; win_dur::Real)
    i = timeindex(t, s.fs)
    k = timeindex(lag, s.fs)
    win_size = timeindex(win_dur, s.fs)
    return ccf(s,i,k;win_size)
end

ccf(s::speech_waveform, i::Int, k::AbstractVector{Int}; win_size::Int) = map(x->ccf(s,i,x;win_size),k)

ccf(s::speech_waveform, t::Real, lags::AbstractVector{<:Real}; win_dur::Real) = map(x->ccf(s,t,x;win_dur),lags)