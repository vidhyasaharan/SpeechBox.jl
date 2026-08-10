#RAPT algorithm parameters from Talkin, D., & Kleijn, W. B. (1995). A robust algorithm for pitch tracking (RAPT). Speech coding and synthesis, 495, 518.
const F0min = 50.0
const F0max = 500.0
const frame_step = 0.01
const nacf_win_size = 0.0075
const CAND_TR = 0.3
const LAG_WT = 0.3
const FREQ_WT = 0.02
const VTRAN_C = 0.005
const VTR_A_C = 0.5
const VTR_S_C = 0.5
const VO_BIAS = 0.0
const DOUBL_C = 0.35
# const A_FACT = 10000.0
const A_FACT = 0.0
const N_CANDS = 20

struct RAPT_candidates
    num_cands::Int
    cands::Vector{Float}
    Φ::Vector{Float}
    index::Int
    time::Float
    fs::Float
end

function RAPT_candidates(num_cands::Int, cands::Vector{<:Real}, Φ::Vector{<:Real}, index::Int, time::Real, fs::Real)
    cands_f = convert(Vector{Float}, cands)
    Φ_f = convert(Vector{Float}, Φ)
    time_f = convert(Float, time)
    fs_f = convert(Float, fs)
    return RAPT_candidates(num_cands, cands_f, Φ_f, index, time_f, fs_f)
end

struct RAPT_timefreq
    pd::timefreq
    candarray::Vector{RAPT_candidates}
end

function RAPT_timefreq(s::speech_waveform; win_dur::Real = nacf_win_size, win_step::Real = frame_step)
    frames = framed_signal(s, win_dur+(1/F0min), win_step)
    pd = amp2db(periodogram(frames))
    candarray = RAPT_pitch_candidates(s; win_dur, win_step)
    return RAPT_timefreq(pd,candarray)
end

function cand_tuples(rtf::RAPT_timefreq)
    ntuples::Int = 0
    lc = RAPT_cands_local_costs(rtf.candarray)
    for cand in rtf.candarray
        ntuples += cand.num_cands
    end
    ptuples = Vector{Tuple{Int, Int}}(undef,ntuples)
    util = Vector{Float}(undef,ntuples)
    # ptuples = Vector{Tuple{Int, Float}}(undef,ntuples)
    i::Int = 1
    for indx ∈ eachindex(rtf.candarray)
        for j in eachindex(rtf.candarray[indx].cands)
            cand = rtf.candarray[indx].cands[j]
            ptuples[i] = (indx, frqindex(cand, rtf.pd.frqs))
            util[i] = lc[j+1,indx]
            # ptuples[i] = (indx, j)
            i+=1
        end
    end
    return ptuples, util
end


#Plot recipe for plotting pitch
@recipe function f(p::RAPT_timefreq; nxticks = 8, nyticks = 8)
    msp = p.pd
    frqs = msp.frqs
    time = msp.time
    xtks = generate_ticks(time,nxticks)
    ytks = generate_ticks(frqs,nyticks)
    ptuples, util = cand_tuples(p)

    @series begin
        seriestype := :heatmap
        xticks := xtks
        yticks := ytks
        # xguide := "Time (sec)"
        # yguide := "Frequency (Hz)"
        # if(~isnothing(msp.title))
        #     title := msp.title
        # end
        msp.components
    end

    @series begin
        seriestype := :scatter
        xticks := xtks
        yticks := ytks
        markersize := util
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


# Estimate pitch candidates from signal at given sample index
function RAPT_pitch_candidates(s::speech_waveform, indx::Int; win_dur::Real = nacf_win_size, ncands::Int = N_CANDS)
    fs = s.fs
    win_size = time2nsamples(win_dur,fs)
    min_lag = time2nsamples(1/F0max,fs)
    max_lag = time2nsamples(1/F0min,fs)
    k = min_lag:max_lag
    cf = nacf(s, indx, k; win_size, nconst = A_FACT)
    inds, mags = findpeaks_sorted(cf)
    thr = CAND_TR*mags[1]
    num_cands = min(sum(mags.>thr), ncands)
    cands = Vector{Float}(undef,num_cands)
    for i ∈ eachindex(cands)
        cands[i] = fs/(min_lag - 1 + inds[i])
    end
    return RAPT_candidates(num_cands,cands,mags[1:num_cands],indx,indx/fs,fs)
end

#Estimate pitch candidates of a signal at every window step (win_step)
function RAPT_pitch_candidates(s::speech_waveform; win_dur::Real = nacf_win_size, win_step::Real = frame_step, ncands::Int = N_CANDS)
    fs = s.fs
    win_size = time2nsamples(win_dur,fs)
    win_shift = time2nsamples(win_step,fs)
    max_lag = time2nsamples(1/F0min,fs)
    nframes = number_signal_frames(s,win_size+max_lag-1,win_shift)
    cand_array = Vector{RAPT_candidates}(undef,nframes)
    for i ∈ eachindex(cand_array)
        indx = (i-1)*win_shift + 1
        cand_array[i] = RAPT_pitch_candidates(s, indx; win_dur, ncands)
    end
    return cand_array
end


function RAPT_maxcands(ca::Vector{RAPT_candidates})
    maxcands::Int = 0
    for array in ca
        maxcands = max(array.num_cands, maxcands)
    end
    return maxcands
end


function RAPT_cands_local_costs(carray::Vector{RAPT_candidates})
    fs = carray[1].fs
    β = LAG_WT/(fs/F0min)

    narrays = length(carray)
    ncands = RAPT_maxcands(carray)
    lcosts = Matrix{Float}(undef,ncands+1,narrays)
    fill!(lcosts,Inf)

    for i ∈ 1:narrays
        lcosts[1,i] = VO_BIAS + maximum(carray[i].Φ)
        for j ∈ 1:carray[i].num_cands
            Cⱼ = carray[i].Φ[j]
            βLⱼ = β/carray[i].cands[j]
            lcosts[j+1,i] = 1 - (Cⱼ*(1-βLⱼ))
        end
    end

    return lcosts
end






function RAPT_nacf(s::speech_waveform; win_dur::Real = nacf_win_size, win_step::Real = frame_step, nconst::Real = A_FACT)
    fs = s.fs
    win_size = time2nsamples(win_dur,fs)
    win_shift = time2nsamples(win_step,fs)
    min_lag = time2nsamples(1/F0max,fs)
    max_lag = time2nsamples(1/F0min,fs)
    nacf(s;min_lag, max_lag, win_size, win_shift, nconst)
end