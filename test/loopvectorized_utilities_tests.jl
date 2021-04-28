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