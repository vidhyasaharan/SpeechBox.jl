
#RAPT algorithm parameters from Talkin, D., & Kleijn, W. B. (1995). A robust algorithm for pitch tracking (RAPT). Speech coding and synthesis, 495, 518.
const F0min = 50.0
const F0max = 500.0
const frame_step = 0.01
const nccf_win_size = 0.0075
const CAND_TR = 0.3
const LAG_WT = 0.3
const FREQ_WT = 0.02
const VTRAN_C = 0.005
const VTR_A_C = 0.5
const VTR_S_C = 0.5
const VO_BIAS = 0.0
const DOUBL_C = 0.35
const A_FACT = 10000.0
const N_CANDS = 20

struct RAPT_candidates
    num_cands::Int
    cands::Vector{Float}
    index::Int
    time::Real
end

struct RAPT_timefreq
    pd::timefreq
    candarray::Vector{RAPT_candidates}
end

function RAPT_timefreq(s::speech_waveform; win_dur::Real = nccf_win_size, win_step::Real = frame_step)
    frames = framed_signal(s, win_dur+(1/F0min), win_step)
    pd = amp2db(periodogram(frames))
    candarray = RAPT_pitch_candidates(s; win_dur, win_step)
    return RAPT_timefreq(pd,candarray)
end

function cand_tuples(rtf::RAPT_timefreq)
    ntuples::Int = 0
    for cand in rtf.candarray
        ntuples += cand.num_cands
    end
    ptuples = Vector{Tuple{Int, Int}}(undef,ntuples)
    i::Int = 1
    for indx ∈ eachindex(rtf.candarray)
        for j in rtf.candarray[indx].cands
            ptuples[i] = (indx, frqindex(j,rtf.pd.frqs))
            i+=1
        end
    end
    return ptuples
end


#Plot recipe for plotting pitch
@recipe function f(p::RAPT_timefreq; nxticks = 8, nyticks = 8)
    msp = p.pd
    frqs = msp.frqs
    time = msp.time
    xtks = generate_ticks(time,nxticks)
    ytks = generate_ticks(frqs,nyticks)
    ptuples = cand_tuples(p)

    @series begin
        seriestype := :heatmap
        xticks := xtks
        yticks := ytks
        msp.components
    end

    @series begin
        seriestype := :scatter
        xticks := xtks
        seriescolor := :green
        legend := false
        xguide := "Time (sec)"
        yguide := "Frequency (Hz)"
        if(~isnothing(msp.title))
            title := msp.title
        end
        ptuples
    end
end



function RAPT_pitch_candidates(s::speech_waveform, i::Int; win_dur::Real = nccf_win_size, ncands::Int = N_CANDS)
    fs = s.fs
    win_size = timeindex(win_dur,fs)
    min_lag = timeindex(1/F0max,fs)
    max_lag = timeindex(1/F0min,fs)
    k = min_lag:max_lag
    cf = nccf(s, i, k; win_size)
    inds, mags = findpeaks_sorted(cf)
    thr = CAND_TR*mags[1]
    num_cands = min(sum(mags.>thr), ncands)
    cands = Vector{Float}(undef,num_cands)
    for i ∈ eachindex(cands)
        cands[i] = fs/(min_lag - 1 + inds[i])
    end
    return RAPT_candidates(num_cands,cands,i,i/fs)
end

function RAPT_pitch_candidates(s::speech_waveform; win_dur::Real = nccf_win_size, win_step::Real = frame_step, ncands::Int = N_CANDS)
    fs = s.fs
    win_size = timeindex(win_dur,fs)
    win_shift = timeindex(win_step,fs)
    max_lag = timeindex(1/F0min,fs)
    nframes = number_signal_frames(s,win_size+max_lag-1,win_shift)
    cand_array = Vector{RAPT_candidates}(undef,nframes)
    for i ∈ eachindex(cand_array)
        indx = (i-1)*win_shift + 1
        cand_array[i] = RAPT_pitch_candidates(s, indx; win_dur, ncands)
    end
    return cand_array
end