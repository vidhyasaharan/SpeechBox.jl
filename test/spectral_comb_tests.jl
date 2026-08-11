@testset "logfreq_spec_comb" begin
    γ = 1.8
    K = 10
    frq_per_octave = 100
    h,z = SpeechBox.generate_logfrq_pitch_comb(;γ, K, frq_per_octave)
    @test typeof(z) == Int
    @test typeof(h) == Vector{Float64}
    pi,pm = SpeechBox.findpeaks(h)
    @test length(pm) == K
    @test pi[1] == z
    noctaves = Int(floor(log2(K)))
    for i=1:noctaves
        @test pi[2^(i-1)] == z + ((i-1)*frq_per_octave)
    end
end

@testset "xcorr_spectral_comb" begin
    fs = 16000
    dur = 2.0
    f = [1000, 5000]
    bw = [200, 500]
    ar = SpeechBox.allpole(f,bw;fs)
    allowed_error = 0.02

    f₀_list = [50 75 100 130 150 170 200 240 280 320 360 410]
    for f₀ in f₀_list
        x = filt(ar, SpeechBox.impulse_train(f₀, dur, fs))
        y, frqs = SpeechBox.xcorr_spectral_comb(x, fs)
        @test typeof(y) == Vector{Float64}
        @test length(y) == length(frqs)
        @test (abs(frqs[argmax(y)]-f₀)/f₀) < allowed_error

        frames = framed_signal(x, fs, 0.10, 0.01)
        y, frqs = SpeechBox.xcorr_spectral_comb(frames)
        @test typeof(y) == Matrix{Float64}
        @test size(y,1) == length(frqs)
        @test size(y,2) == frames.num_signal_frames
        for i=1:size(y,2)
            @test (abs(frqs[argmax(y[:,i])]-f₀)/f₀) < allowed_error
        end
    end


end