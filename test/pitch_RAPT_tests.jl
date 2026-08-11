
@testset "RAPT_candidates" begin
    fs = 16000
    ncands = 10
    rcands = rand(Float64,ncands)
    mags = rand(Float64,ncands)
    time = 0.04
    indx = Int(round(time*fs))
    cands = SpeechBox.RAPT_candidates(length(rcands), rcands, mags, indx, time, fs)
    @test cands isa SpeechBox.RAPT_candidates{Float64}
    @test typeof(cands.num_cands) == Int
    @test typeof(cands.cands) == Vector{Float64}
    @test typeof(cands.Φ) == Vector{Float64}
    @test typeof(cands.index) == Int
    @test typeof(cands.time)<:Real
    @test length(cands.cands) == cands.num_cands
    @test length(cands.Φ) == cands.num_cands
end


@testset "RAPT_pitch_candidates_index" begin
    indx = Int(round((10*SpeechBox.frame_step*sig.fs) + 1))
    cands = SpeechBox.RAPT_pitch_candidates(sig, indx)
    @test cands isa SpeechBox.RAPT_candidates{Float64}
    @test cands.num_cands == length(cands.cands)
    @test cands.index == indx
    @test cands.time == indx/sig.fs
end

@testset "RAPT_pitch_candidates" begin
    cand_array = SpeechBox.RAPT_pitch_candidates(sig)
    @test typeof(cand_array) == Vector{SpeechBox.RAPT_candidates{Float64}}
end

@testset "RAPT_maxcands" begin
    cand_array = SpeechBox.RAPT_pitch_candidates(sig)
    mxcands = SpeechBox.RAPT_maxcands(cand_array)
    @test typeof(mxcands) == Int
    @test mxcands == SpeechBox.N_CANDS

    cand_array = SpeechBox.RAPT_pitch_candidates(sig; ncands = 10)
    mxcands = SpeechBox.RAPT_maxcands(cand_array)
    @test mxcands == 10
end