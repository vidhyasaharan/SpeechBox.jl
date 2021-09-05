using Test
using SpeechBox
using LinearAlgebra
using WAV
using DSP: filt
using SpeechBox: Float
using Statistics

testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"))
srcdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../src/"))

x, fs = wavread(joinpath(testdir,"King.wav"))
x = x[:]
signal = speech_waveform(x,fs)
