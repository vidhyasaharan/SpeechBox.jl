@testset "dotavx" begin
    ndim = 10
    a = rand(Float,ndim)
    b = rand(Float,ndim)
    ca = rand(Complex{Float},ndim)
    cb = rand(Complex{Float},ndim)
    @test a⋅b ≈ SpeechBox.dotavx(a,b)
    @test a⋅cb ≈ SpeechBox.dotavx(a,cb)
    @test ca⋅b ≈ SpeechBox.dotavx(ca,b)
    @test ca⋅cb ≈ SpeechBox.dotavx(ca,cb)
end

@testset "xcorr" begin
    len = 1000
    hlen = 100
    hstart = 101
    x = randn(Float,len)
    h = x[hstart:hstart+hlen-1]
    z = 10
    y = SpeechBox.xcorr(x,h,z)
    @test typeof(y) == Vector{Float}
    @test length(y) == length(x)
    @test argmax(y) == hstart + z - 1
    @test maximum(y) ≈ SpeechBox.dotavx(h)
end