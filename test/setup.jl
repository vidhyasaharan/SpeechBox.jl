using Test
using SpeechBox
using LinearAlgebra
using WAV
using DSP: filt
using Statistics
using Random

#Fix the RNG seed so the stochastic tests (k-means/GMM initialisation and sampling) are
#deterministic - each @testset re-seeds from this state. white_noise/ar_process draw from a
#freshly entropy-seeded generator by design and stay nondeterministic; the tests that use
#them (lpc) only assert loose statistical tolerances.
Random.seed!(2026)

const Float = Float64 #Concrete type used by the test fixtures

testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"))
srcdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../src/"))

x, fs = wavread(joinpath(testdir,"King.wav"))
x = x[:]
signal = speech_waveform(x,fs)
