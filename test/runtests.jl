

using Test
using SpeechBox
using LinearAlgebra
# using LibSndFile
# using FileIO
# using SampledSignals
using WAV
using DSP: filt
using SpeechBox: Float



testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"))
srcdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../src/"))

x, fs = wavread(joinpath(testdir,"King.wav"))
x = x[:]
signal = speech_waveform(x,fs)

@testset verbose = true "SpeechBox" begin
    @testset "utilities" begin 
        include("internal_utilities_tests.jl")
        include("loopvectorized_utilities_tests.jl")
    end
    @testset "correlation" begin
        include("correlation_tests.jl")
    end
    @testset "structs" begin
        include("signal_objects_tests.jl")
    end
    @testset "spectrum" begin 
        include("spectralanalyses_tests.jl")
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
end