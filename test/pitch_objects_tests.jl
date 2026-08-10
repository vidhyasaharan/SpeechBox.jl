@testset "pitch_timefreq" begin
    #Half a second of the test recording keeps this fast
    s = speech_waveform(signal.x[1:8000], signal.fs)
    frames = framed_signal(s, 0.02, 0.01)
    p = pitch(spectral_comb(), frames)
    ptf = SpeechBox.pitch_timefreq(p, frames)

    @test ptf isa SpeechBox.pitch_timefreq{Float,Float}
    @test length(ptf.pitch) == length(p)
    @test length(ptf.pindx) == length(p)
    @test ptf.msp isa timefreq

    #Voiced frames map to the nearest spectrogram row, unvoiced frames stay NaN
    row_spacing = ptf.msp.frqs[2] - ptf.msp.frqs[1]
    for i ∈ eachindex(p)
        if isnan(p[i])
            @test isnan(ptf.pindx[i])
        else
            @test 1 <= ptf.pindx[i] <= length(ptf.msp.frqs)
            @test abs(ptf.msp.frqs[Int(ptf.pindx[i])] - p[i]) <= row_spacing
        end
    end

    #A pitch contour of the wrong length is rejected
    @test SpeechBox.pitch_timefreq(p[1:end-1], frames) === nothing
end
