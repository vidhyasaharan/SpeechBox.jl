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
