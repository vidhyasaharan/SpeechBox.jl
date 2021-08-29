@testset "running_mean" begin
    ndim = 2
    npts = 100000
    x = randn(Float,ndim,npts)
    @test mean(x,dims=2) ≈ SpeechBox.running_mean(x)
end

@testset "running_meanvar" begin
    ndim = 2
    npts = 100000
    x = randn(Float,ndim,npts)
    m,v = SpeechBox.running_meanvar(x)
    @test mean(x,dims=2) ≈ m
    @test var(x,dims=2) ≈ v
end