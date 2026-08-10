@testset "window" begin
    flen = 5
    fmid = Int(ceil(flen/2))

    win = SpeechBox.window(flen;wtype="rect")
    @test typeof(win)<:Array{Float,1}
    @test length(win) == flen
    for i=1:flen
        @test win[i] == 1.0
    end


    win = SpeechBox.window(flen;wtype="hanning")
    @test typeof(win)<:Array{Float,1}
    @test length(win) == flen
    α₀ = 0.5
    α₁ = 1 - α₀
    N = flen-1
    wfun(n) = α₀ - α₁*cos(2π*n/N)
    for i=1:flen
        @test win[i] ≈ wfun(i-1)
    end


    win = SpeechBox.window(flen;wtype="hamming")
    @test typeof(win)<:Array{Float,1}
    @test length(win) == flen
    α₀ = 0.54
    α₁ = 1 - α₀
    N = flen-1
    wfun1(n) = α₀ - α₁*cos(2π*n/N)
    for i=1:flen
        @test win[i] ≈ wfun1(i-1)
    end


    win1 = SpeechBox.window(flen)
    win2 = SpeechBox.window(flen;wtype="hanning")
    @test win1 == win2
end


@testset "logfreq_array" begin
    frqs = SpeechBox.logfreq_array(;fmin = 2, fmax = 8192, frq_per_octave = 100)
    lfrqs = log2.(frqs)
    @test lfrqs[1] == 1.0
    @test lfrqs[end] == 13.0
    @test length(frqs) == 1201
    for i=1:length(frqs)
        @test frqs[i] ≈ exp2(0.99 + 0.01*i)
    end
end

@testset "linfreq_array" begin
    frqs = SpeechBox.linfreq_array(;fmin = 0, fmax = 100, nfrqs = 11)
    @test frqs[1] == 0.0
    @test frqs[end] == 100.0
    for i=2:length(frqs)
        @test frqs[i] == (i-1)*10
    end
end


@testset "cexp" begin
    fs = 8000
    N = 16000
    cs(fr,i) = cos(2π*(fr/fs)*i)
    sn(fr,i) = sin(2π*(fr/fs)*i)

    frqs = [100,700,3401]

    for f in frqs
        f = 100
        ce = sqrt(N)*SpeechBox.cexp(f,fs,N)
        @test ce⋅ce ≈ N
        for i=1:N
            @test isapprox(real(ce[i]),cs(f,i);atol = 1e-10)
            @test isapprox(imag(ce[i]),-sn(f,i);atol = 1e-10)
        end
    end
end

@testset "cexp_proj_matrix" begin
    fs = 16000
    frqs = [10, 50, 100, 500, 1000, 5000]
    N = 1600
    proj = SpeechBox.cexp_proj_matrix(frqs,fs,N)
    @test typeof(proj) <: Array{Complex{Float},2}
    @test size(proj,1) == length(frqs)
    @test size(proj,2) == N
    for i = 1:length(frqs)
        ce = SpeechBox.cexp(frqs[i],fs,N)
        @test proj[i,:] == ce
    end
end


@testset "Δ" begin
    x = [1,1,2,3,4,5]
    y = SpeechBox.Δ(x)
    @test typeof(y) <: AbstractVector
    @test length(y) == length(x) - 1
    @test y[1] == 0
    for i=2:length(y)
        @test y[i] == 1
    end
end

@testset "findpeaks" begin
    x = [3,0,4,2,1,2,3,2,1]

    ind,mag = SpeechBox.findpeaks(x;min_dist = 1)
    @test typeof(ind) == Vector{Int}
    @test typeof(mag) == Vector{Int}
    @test length(ind) == length(mag)
    @test length(ind) == 3
    for i in eachindex(ind)
        @test x[ind[i]] == mag[i]
    end


end


@testset "remove_nearest_peak!" begin
    x = [3,0,4,2,1,2,3,2,1]
    ind,mag = SpeechBox.findpeaks(x;min_dist = 1)
    SpeechBox.remove_nearest_peak!(ind,mag)
    @test typeof(ind) == Vector{Int}
    @test typeof(mag) == Vector{Int}
    @test length(ind) == 2
    @test length(mag) == 2
    @test ind[1] == 3
    @test ind[2] == 7
    @test mag[1] == 4
    @test mag[2] == 3
end


@testset "findpeaks-dist" begin
    x = [3,0,4,2,1,2,3,2,1]

    ind,mag = SpeechBox.findpeaks(x;min_dist = 3)
    @test length(ind) == 2
    @test ind[1] == 3
    @test ind[2] == 7
end

@testset "findpeaks_sorted" begin
    x = [3,0,4,2,1,2,3,2,1]

    ind,mag = SpeechBox.findpeaks_sorted(x;min_dist = 1)
    @test typeof(ind) == Vector{Int}
    @test typeof(mag) == Vector{Int}
    @test length(ind) == length(mag)
    @test length(ind) == 3
    @test mag[1] == 4
    @test mag[2] == 3
    @test mag[3] == 3
    @test ind[1] == 3
    @test ind[2] == 1
    @test ind[3] == 7
end


@testset "resample" begin

    fs = 8000
    n = 1:2*fs
    frq = 1000
    xx = cos.(2*pi*(frq/fs)*n)

    s = SpeechBox.speech_waveform(xx,fs)
    rs = SpeechBox.resample(s,2*fs)
    @test rs isa SpeechBox.speech_waveform{Float}
    @test abs(length(rs.x) - length(s.x)*2) <= 1
    @test rs.fs == s.fs*2

    mag = magspec(s)
    rmag = magspec(rs)

    ind = argmax(mag.components)
    rind = argmax(rmag.components)

    @test abs(mag.frqs[ind] - rmag.frqs[rind]) < abs(mag.frqs[ind] - rmag.frqs[rind-1])
    @test abs(mag.frqs[ind] - rmag.frqs[rind]) < abs(mag.frqs[ind] - rmag.frqs[rind+1])

end

@testset "findclosest" begin
    data = 0.0:0.1:10.0
    rin = 7
    x = (rin-1)*0.1
    cin = SpeechBox.findclosest(x,data)
    @test cin == rin
end


@testset "frqindex" begin
    frqs = 0.0:10.0:1000.0
    f = 30.0
    fin = SpeechBox.frqindex(f,frqs)
    @test typeof(fin) == Int
    @test fin == 4

    fv = [21.0, 41.5, 59.0]
    fvin = SpeechBox.frqindex(fv,frqs)
    @test typeof(fvin) == Vector{Int}
    @test length(fvin) == length(fv)
    @test fvin[1] == 3
    @test fvin[2] == 5
    @test fvin[3] == 7
end


@testset "element_op" begin
    e = SpeechBox.element_op_spectrum("Base.log10")
    @test typeof(e) == Expr
    e = SpeechBox.element_op_timefreq("Base.log10")
    @test typeof(e) == Expr
    frames = SpeechBox.framed_signal(x,fs,0.09,0.01)
    frame = SpeechBox.extract_frame(frames,11)
    mag = SpeechBox.magspec(frame,fs)
    msp = SpeechBox.specgram(frames)
    lmag = log(mag)
    l10mag = log10(mag)
    lmsp = log(msp)
    l10msp = log10(msp)
    @test lmag.components == log.(mag.components)
    @test l10mag.components == log10.(mag.components)
    @test lmsp.components == log.(msp.components)
    @test l10msp.components == log10.(msp.components)
end
