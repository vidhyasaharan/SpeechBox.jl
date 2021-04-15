
using Test
using SpeechBox
using LinearAlgebra
# using LibSndFile
# using FileIO
# using SampledSignals
using WAV
using Random

testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"))
srcdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../src/"))

x, fs = wavread(joinpath(testdir,"King.wav"))
x = x[:]
signal = speech_waveform(x,fs)

include("internal_utilities_tests.jl")
include("signal_objects_tests.jl")
include("spectralanalyses_tests.jl")
include("mfcc_tests.jl")
include("vad_tests.jl")
include("lpc_tests.jl")
