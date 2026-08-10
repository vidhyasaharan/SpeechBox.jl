@testset "Δ" begin
    x = [1,1,2,3,4,5]
    y = SpeechBox.Δ(x)
    @test typeof(y) <: AbstractVector
    @test length(y) == length(x) - 1
    @test y[1] == 0
    for i=2:length(y)
        @test y[i] == 1
    end
end

@testset "findpeaks" begin
    x = [3,0,4,2,1,2,3,2,1]

    ind,mag = SpeechBox.findpeaks(x;min_dist = 1)
    @test typeof(ind) == Vector{Int}
    @test typeof(mag) == Vector{Int}
    @test length(ind) == length(mag)
    @test length(ind) == 3
    for i in eachindex(ind)
        @test x[ind[i]] == mag[i]
    end
end


@testset "remove_nearest_peak!" begin
    x = [3,0,4,2,1,2,3,2,1]
    ind,mag = SpeechBox.findpeaks(x;min_dist = 1)
    SpeechBox.remove_nearest_peak!(ind,mag)
    @test typeof(ind) == Vector{Int}
    @test typeof(mag) == Vector{Int}
    @test length(ind) == 2
    @test length(mag) == 2
    @test ind[1] == 3
    @test ind[2] == 7
    @test mag[1] == 4
    @test mag[2] == 3
end


@testset "findpeaks-dist" begin
    x = [3,0,4,2,1,2,3,2,1]

    ind,mag = SpeechBox.findpeaks(x;min_dist = 3)
    @test length(ind) == 2
    @test ind[1] == 3
    @test ind[2] == 7
end

@testset "findpeaks_sorted" begin
    x = [3,0,4,2,1,2,3,2,1]

    ind,mag = SpeechBox.findpeaks_sorted(x;min_dist = 1)
    @test typeof(ind) == Vector{Int}
    @test typeof(mag) == Vector{Int}
    @test length(ind) == length(mag)
    @test length(ind) == 3
    @test mag[1] == 4
    @test mag[2] == 3
    @test mag[3] == 3
    @test ind[1] == 3
    @test ind[2] == 1
    @test ind[3] == 7
end
