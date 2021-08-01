
@testset "RAPT_candidates" begin
    fs = 16000
    ncands = 10
    rcands = rand(Float,ncands)
    mags = rand(Float,ncands)
    time = 0.04
    indx = Int(round(time*fs))
    cands = SpeechBox.RAPT_candidates(length(rcands), rcands, mags, indx, time, fs)
    @test typeof(cands) == SpeechBox.RAPT_candidates
    @test typeof(cands.num_cands) == Int
    @test typeof(cands.cands) == Vector{Float}
    @test typeof(cands.Φ) == Vector{Float}
    @test typeof(cands.index) == Int
    @test typeof(cands.time)<:Real
    @test length(cands.cands) == cands.num_cands
    @test length(cands.Φ) == cands.num_cands
end


@testset "RAPT_pitch_candidates_index" begin
    indx = Int(round((10*SpeechBox.frame_step*signal.fs) + 1))
    cands = SpeechBox.RAPT_pitch_candidates(signal, indx)
    @test typeof(cands) == SpeechBox.RAPT_candidates
    @test cands.num_cands == length(cands.cands)
    @test cands.index == indx
    @test cands.time == indx/signal.fs
end

@testset "RAPT_pitch_candidates" begin
    cand_array = SpeechBox.RAPT_pitch_candidates(signal)
    @test typeof(cand_array) == Vector{SpeechBox.RAPT_candidates}
end

@testset "RAPT_maxcands" begin
    cand_array = SpeechBox.RAPT_pitch_candidates(signal)
    mxcands = SpeechBox.RAPT_maxcands(cand_array)
    @test typeof(mxcands) == Int
    @test mxcands == SpeechBox.N_CANDS

    cand_array = SpeechBox.RAPT_pitch_candidates(signal; ncands = 10)
    mxcands = SpeechBox.RAPT_maxcands(cand_array)
    @test mxcands == 10
end