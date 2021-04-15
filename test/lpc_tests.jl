@testset "lpc_order" begin
    fs = 8000
    N = SpeechBox.lpc_order(fs)
    @test N == 10
    @test SpeechBox.lpc_order(16000.25) == 18
end


@testset "lpc" begin
    c = 0.8
    d = 3
    N = 5
    a = zeros(d+1)
    a[1] = 1.0
    a[d+1] = c

    x = SpeechBox.ar_process(a, 10000)

    coeff = lpc(x,N)

    @test length(coeff) == N+1
    @test coeff[1] == 1.0
    
    @test coeff[d+1] < 1.1*c
    @test coeff[d+1] > 0.9*c
    for i=2:length(coeff)
        if i≠(d+1)
            @test abs(coeff[i]) < 0.1*c
        end
    end
end


@testset "lpc_freqz" begin
    fs = 16000
    nsam = 10000
    ar,f,bw = SpeechBox.rand_allpole(fs, 1)
    x = filt(ar, SpeechBox.white_noise(nsam))
    
    nfrqs = 10
    frqs = [f; rand(nfrqs)*(fs/2)]
    N = 4
    h = SpeechBox.lpc_freqz(x,fs,N;frqs)
    for i=2:length(h)
        @test abs(h[1]) > abs(h[i])
    end
end