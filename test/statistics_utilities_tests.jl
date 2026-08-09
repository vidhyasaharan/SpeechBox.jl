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

@testset "running_meancov" begin
    ndim = 20
    npts = 10000
    x = randn(Float64,ndim,npts)
    m,C = SpeechBox.running_meancov(x)
    m̂ = mean(x,dims=2)
    Ĉ = ((x.-m̂)*(x.-m̂)')/(size(x,2)-1)
    @test m̂ ≈ m
    @test Ĉ ≈ C
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

@testset "distance computations" begin
    ndim = 10
    a = rand(Float,ndim)
    b = rand(Float,ndim)
    @test (a-b)⋅(a-b) ≈ SpeechBox.sqL2dist(a,b)
    @test sqrt(sum(abs2,a-b)) ≈ SpeechBox.L2dist(a,b)

    A = randn(Float,ndim,ndim)
    x = randn(Float,ndim)
    y = randn(Float,ndim)
    d = (x-y)'*A*(x-y)
    @test SpeechBox.sqmahal(x,y,A) ≈ d
end

@testset "Distance" begin
    a = randn(Float64,10)
    b = randn(Float64,10)
    c1 = [a b a b]
    c2 = [a b]

    @test SpeechBox.SqL2 <: SpeechBox.Distance
    @test SpeechBox.L2 <: SpeechBox.Distance

    @test SpeechBox.dist(SpeechBox.SqL2(),a,b) ≈ SpeechBox.sqL2dist(a,b)
    @test SpeechBox.dist(SpeechBox.L2(),a,b) ≈ SpeechBox.L2dist(a,b)

    distances = (SpeechBox.SqL2(), SpeechBox.L2())
    for dm ∈ distances
        d = SpeechBox.pairwise(dm,a,c1)
        @test length(d) == size(c1,2)
        @test d[1] == 0.0
        @test d[2] == SpeechBox.dist(dm, a, b)

        dd = SpeechBox.pairwise(dm, c2, c1)
        @test size(dd,1) == size(c2,2)
        @test size(dd,2) == size(c1,2)
        @test dd[1,1] == 0.0
        @test dd[2,2] == 0.0
        @test dd[1,2] == SpeechBox.dist(dm, a, b)
        @test dd[2,1] == SpeechBox.dist(dm, b, a)
    end

end