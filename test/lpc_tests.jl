@testset "lpc_order" begin
    fs = 8000
    N = SpeechBox.lpc_order(fs)
    @test N == 10
    @test SpeechBox.lpc_order(16000.25) == 18
end


