
using Test
using SpeechBox
# using LibSndFile
# using FileIO
# using SampledSignals
using WAV

testdir = normpath(joinpath(dirname(@__FILE__),"../test/"))
srcdir = normpath(joinpath(dirname(@__FILE__),"../src/"))

x, fs = wavread(joinpath(testdir,"King.wav"))

include("speechframing_tests.jl")
include("spectralanalyses_tests.jl")
include("mfcc_tests.jl")
include("vad_tests.jl")
