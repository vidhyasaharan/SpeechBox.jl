
include("setup.jl")

@testset verbose = true "SpeechBox" begin
    @testset "utilities" begin
        include("internal_utilities_tests.jl")
        include("statistics_utilities_tests.jl")
    end
    @testset "correlation" begin
        include("correlation_tests.jl")
    end
    @testset "structs" begin
        include("pitch_objects_tests.jl")
    end
    @testset "spectrum" begin
        include("lpc_tests.jl")
    end
    @testset "vad" begin
        include("vad_tests.jl")
    end
    @testset "MFCC" begin
        include("mfcc_tests.jl")
    end
    @testset "spec. comb" begin
        include("spectral_comb_tests.jl")
    end
    @testset "RAPT" begin
        include("pitch_RAPT_tests.jl")
    end
    @testset "k-means" begin
        include("kmeans_tests.jl")
    end
    @testset "GMM" begin
        include("GMM_tests.jl")
    end
    @testset "Float32" begin
        include("float32_tests.jl")
    end
end
