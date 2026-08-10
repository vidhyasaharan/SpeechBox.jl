@testset "magspec" begin
    t = 0:0.01:0.99
    frq = 25
    fs = 100.0
    xx = cos.(2*pi*frq*t)
    spec = magspec(xx,fs;wtype = "rect")
    mspec1, frqs1 = magspec(comp(), xx, fs; wtype = "rect")
    mspec = spec.components
    @test mspec1 == mspec1
    @test frqs1 == spec.frqs
    mmag,mfrq = findmax(mspec)
    @test length(mspec) == 51
    @test mfrq == frq + 1
    @test round(mmag) == 50.0
end

@testset "periodogram" begin
    fs = 8000
    n = 1:160
    frq = 1024
    xx = cos.(2*pi*(frq/fs)*n)

    tfrqs = [100, 500, 705, 1024, 1800, 2100, 3401, 3700]
    spec = periodogram(xx,fs,tfrqs)
    @test typeof(spec) == SpeechBox.spectrum{Float}
    pgram  = spec.components
    @test pgram == periodogram(comp(), xx,fs,tfrqs)
    @test length(pgram) == length(tfrqs)
    mmag, mfindx = findmax(pgram)
    @test tfrqs[mfindx] == frq

    spec = periodogram(xx,fs;fmin=8,fmax=4000)
    @test typeof(spec) == SpeechBox.spectrum{Float}
    pgram = spec.components
    @test pgram == periodogram(comp(),xx,fs,fmin=8,fmax=4000)
    mmag, mfindx = findmax(pgram)
    frqs = SpeechBox.logfreq_array(;fmin=8,fmax=4000)
    @test length(pgram) == length(frqs)
    @test frqs[mfindx] == frq
end


@testset "specgram" begin
    tf = specgram(x, fs;frame_dur = 0.03,frame_shift_dur=0.01, wtype = "hamming")
    msp = tf.components
    @test msp == specgram(comp(), x, fs;frame_dur = 0.03,frame_shift_dur=0.01, wtype = "hamming")
    @test typeof(msp) == Matrix{Float}
    @test size(msp,1) > 0
    @test size(msp,2) > 0
    @test length(tf.frqs) == size(msp,1) #one frequency index per spectrogram row
    @test tf.frqs[1] == 0.0
    @test tf.frqs[end] <= fs/2
    @test maximum(isa.(msp,Complex))==false
    @test minimum(msp) >= eps()

    t = 0:0.01:0.99
    frq1 = 25
    frq2 = 30
    xx = [cos.(2*pi*frq1*t);cos.(2*pi*frq2*t)]
    tf = specgram(xx,100.0,frame_dur = 1.0, frame_shift_dur = 1.0)
    msp = tf.components
    @test size(msp,1) > 0
    @test size(msp,2) > 0
    @test maximum(isa.(msp,Complex))==false
    @test minimum(msp) >= eps()
    mmag1,mfrq1 = findmax(msp[:,1])
    mmag2,mfrq2 = findmax(msp[:,2])
    @test mfrq1 == frq1 + 1
    @test mfrq2 == frq2 + 1
    @test length(tf.frqs) == size(msp,1)
    @test tf.frqs[mfrq1] ≈ frq1 #frequency axis maps spectrogram rows to Hz
    @test tf.frqs[mfrq2] ≈ frq2

end
