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