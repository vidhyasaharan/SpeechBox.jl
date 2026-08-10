#Verify that the processing chain is generic over the floating point precision of the input signal
@testset "objects" begin
    x32 = randn(Float32, 16000)
    s32 = speech_waveform(x32, 16000)
    @test s32 isa speech_waveform{Float32}
    @test typeof(s32.fs) == Float32

    frames32 = framed_signal(s32)
    @test frames32 isa framed_signal{Float32}
    @test eltype(extract_frame(frames32,1)) == Float32
    @test eltype(frame_energy(frames32)) == Float32
end

@testset "framing and windows" begin
    x32 = randn(Float32, 1600)
    @test eltype(SpeechBox.window(Float32, 32)) == Float32
    @test SpeechBox.window(Float32, 32) ≈ SpeechBox.window(32)
    @test eltype(enframe(x32, 160, 80)) == Float32
    @test eltype(preemphasis(x32)) == Float32
end

@testset "spectral analyses" begin
    x32 = randn(Float32, 16000)
    s32 = speech_waveform(x32, 16000)
    frames32 = framed_signal(s32)

    @test eltype(dft(randn(Float32,256))) == Complex{Float32}

    ms = magspec(s32)
    @test ms isa SpeechBox.spectrum{Float32,Float32}
    @test eltype(ms.frqs) == Float32

    tf = specgram(frames32)
    @test tf isa SpeechBox.timefreq{Float32,Float32}
    @test length(tf.frqs) == size(tf.components,1)
    @test eltype(tf.time) == Float32

    pg = periodogram(SpeechBox.comp(), frames32)
    @test eltype(pg) == Float32

    #Float32 and Float64 paths should agree numerically
    x64 = convert(Vector{Float64}, x32)
    tf64 = specgram(framed_signal(x64, 16000.0))
    @test isapprox(tf.components, tf64.components; rtol = 1e-3)
end

@testset "speech analyses" begin
    x32 = randn(Float32, 16000)
    s32 = speech_waveform(x32, 16000)
    frames32 = framed_signal(s32)

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
