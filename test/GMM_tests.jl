@testset "normaliseWeights" begin
    nmix = 10
    w = rand(Float,nmix)
    ŵ = copy(w)
    SpeechBox.normaliseWeights!(w)
    @test sum(w) ≈ one(Float)
    for i ∈ eachindex(w,ŵ)
        @test w[i]/sum(w) ≈ ŵ[i]/sum(ŵ)
    end
end

# @testset "isconsistentX" begin
#     nmix = 5
#     ndim = 3
#     w = rand(nmix)
# end