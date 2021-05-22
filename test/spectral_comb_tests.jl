@testset "logfreq_spec_comb" begin
    γ = 1.8
    K = 5
    frq_per_octave = 1000
    h,z = SpeechBox.generate_logfrq_pitch_comb(;γ, K, frq_per_octave)
    frqs = SpeechBox.logfreq_array(fmin = 0.5, fmax = K + 0.5, frq_per_octave = frq_per_octave)
    pi,pm = SpeechBox.findpeaks(h)
    @test length(pm) == K
    @test pi[1] == z
    # for i=1:K
    #     @test frqs[pi[i]] ≈ i
    # end
end