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


@testset "lpc_response" begin
    fs = 16000
    nsam = 2*fs
    f = [1000, 5000]
    bw = [200, 200]
    ar = SpeechBox.allpole(f,bw;fs)
    x = filt(ar, SpeechBox.white_noise(nsam))

    freqs = collect(0:100:fs/2)
    lsp = lpc_response(x,fs,SpeechBox.lpc_order(fs); frqs = freqs)

    i,m = SpeechBox.findpeaks(lsp.components)
    
    @test lsp isa SpeechBox.spectrum{Float,Float}
    @test lsp.components == lpc_response(comp(),x,fs,SpeechBox.lpc_order(fs); frqs = freqs)
    @test length(lsp.frqs) == length(freqs)
    @test lsp.frqs[i[1]] == f[1]
    @test lsp.frqs[i[2]] == f[2]

    frames = framed_signal(signal,0.09,0.01)
    lspec = lpc_response(frames, SpeechBox.lpc_order(signal.fs); frqs = freqs)

    @test lspec isa SpeechBox.timefreq{Float,Float}
    @test lspec.components == lpc_response(comp(), frames, SpeechBox.lpc_order(signal.fs); frqs = freqs)
    @test size(lspec.components,1) == length(freqs)
    @test size(lspec.components,2) == frames.num_signal_frames
end
