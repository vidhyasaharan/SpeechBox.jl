@testset "xcorr" begin
    len = 1000
    hlen = 100
    hstart = 101
    x = randn(Float,len)
    h = x[hstart:hstart+hlen-1]
    z = 10
    y = SpeechBox.xcorr(x,h,z)
    y1 = Vector{Float}(undef,length(x))
    SpeechBox.xcorr!(y1, x, h, z)
    @test typeof(y) == Vector{Float}
    @test length(y) == length(x)
    @test y1 == y
    @test argmax(y) == hstart + z - 1
    @test maximum(y) ≈ sum(abs2,h)
end

@testset "acorr" begin
    len = 1000
    p = 100
    x = randn(Float, len)
    rxx = Vector{Float}(undef,p)
    SpeechBox.acorr!(rxx,x)
    rx = SpeechBox.acorr(x,p)
    @test typeof(rx) == Vector{Float}
    @test length(rx) == p
    @test rxx == rx
    @test argmax(rx) == 1
    @test maximum(rx) ≈ sum(abs2,x)/len
end