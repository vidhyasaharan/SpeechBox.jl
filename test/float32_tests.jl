#Verify the speech analyses stay in the precision of the input signal. (The time-frequency
#core has its own Float32 suite in TimeFrequencyAnalysis.jl.)
@testset "waveform alias" begin
    @test speech_waveform === waveform #backwards-compatible name for the TFA container
    s32 = speech_waveform(randn(Float32,8000), 8000)
    @test s32 isa speech_waveform{Float32}
end

@testset "speech analyses" begin
    x32 = randn(Float32, 16000)
    s32 = speech_waveform(x32, 16000)
    frames32 = framed_signal(s32)

    @test eltype(preemphasis(x32)) == Float32
    @test eltype(melfcc(frames32)) == Float32
    @test length(vad(frames32)) == frames32.num_frames
    @test eltype(lpc(randn(Float32, 4000), 10)) == Float32
    @test eltype(SpeechBox.xcorr(randn(Float32,100), randn(Float32,10))) == Float32
    @test eltype(acorr(randn(Float32,1000), 10)) == Float32
    @test eltype(pitch(spectral_comb(), frames32)) == Float32
end

@testset "clustering and GMM" begin
    data32 = 0.25f0 .* randn(Float32, 2, 400) .+ repeat(Float32[1 1 -1 -1; 1 -1 1 -1], 1, 100)
    c32 = kmeans(data32, 4)
    @test eltype(c32) == Float32

    G32, LL = SpeechBox.trainML(data32, 2, 3)
    @test G32 isa SpeechBox.GMM{Float32}
    @test eltype(LL) == Float32
end
