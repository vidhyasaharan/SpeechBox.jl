@testset "running_mean" begin
    ndim = 20
    npts = 10000
    x = randn(Float,ndim,npts)
    @test mean(x,dims=2) ≈ SpeechBox.running_mean(x)
end

@testset "running_meanvar" begin
    ndim = 20
    npts = 10000
    x = randn(Float,ndim,npts)
    m,v = SpeechBox.running_meanvar(x)
    @test mean(x,dims=2) ≈ m
    @test var(x,dims=2) ≈ v
end

@testset "update_running_mean!" begin
    ndim = 20
    npts = 10000
    x = randn(Float,ndim,npts)
    m = zeros(Float,ndim)
    for i ∈ axes(x,2)
        SpeechBox.update_running_mean!(m,x[:,i],i)
    end
    @test mean(x,dims=2) ≈ m
end