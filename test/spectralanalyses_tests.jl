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
    @test length(pgram) == length(tfrqs)
    mmag, mfindx = findmax(pgram)
    @test tfrqs[mfindx] == frq

    spec = periodogram(xx,fs;fmin=8,fmax=4000)
    @test typeof(spec) == SpeechBox.spectrum{Float}
    pgram = spec.components
    mmag, mfindx = findmax(pgram)
    frqs = SpeechBox.logfreq_array(;fmin=8,fmax=4000)
    @test length(pgram) == length(frqs)
    @test frqs[mfindx] == frq
end


@testset "specgram" begin
    tf = specgram(x, fs;win_dur = 0.03,win_shift=0.01, wtype = "hamming")
    msp = tf.components
    @test typeof(msp) <: Array{<:AbstractFloat,2}
    @test size(msp,1) > 0
    @test size(msp,2) > 0
    @test maximum(isa.(msp,Complex))==false
    @test minimum(msp) >= eps()

    t = 0:0.01:0.99
    frq1 = 25
    frq2 = 30
    xx = [cos.(2*pi*frq1*t);cos.(2*pi*frq2*t)]
    tf = specgram(xx,100.0,win_dur = 1.0, win_shift = 1.0)
    msp = tf.components
    @test size(msp,1) > 0
    @test size(msp,2) > 0
    @test maximum(isa.(msp,Complex))==false
    @test minimum(msp) >= eps()
    mmag1,mfrq1 = findmax(msp[:,1])
    mmag2,mfrq2 = findmax(msp[:,2])
    @test mfrq1 == frq1 + 1
    @test mfrq2 == frq2 + 1

end
