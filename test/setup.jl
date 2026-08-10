using Test
using SpeechBox
using LinearAlgebra
using WAV
using DSP: filt
using Statistics

const Float = Float64 #Concrete type used by the test fixtures

testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"))
srcdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../src/"))

x, fs = wavread(joinpath(testdir,"King.wav"))
x = x[:]
signal = speech_waveform(x,fs)
