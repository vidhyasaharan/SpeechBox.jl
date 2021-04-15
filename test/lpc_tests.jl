@testset "lpc_order" begin
    fs = 8000
    N = SpeechBox.lpc_order(fs)
    @test N == 10
    @test SpeechBox.lpc_order(16000.25) == 18
end


@testset "lpc" begin
    rng = MersenneTwister(100)
    x = randn(rng,Float64,(10000))
    a = 0.8
    d = 3
    N = 5
    for i=d+1:length(x)
        x[i] += a*x[i-d]
    end
    coeff = lpc(x,N)
    @test length(coeff) == N+1
    @test coeff[1] == 1.0
    
    @test coeff[d+1] > -1.1*a
    @test coeff[d+1] < -0.9*a
    for i=2:length(coeff)
        if i≠(d+1)
            @test abs(coeff[i]) < 0.1*a
        end
    end
end