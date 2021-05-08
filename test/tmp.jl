using SpeechBox
using WAV
using Plots
using BenchmarkTools

testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"))
x, fs = wavread(joinpath(testdir,"King.wav"))

frames = framed_signal(x,fs,0.05,0.01)
frame = extract_frame(frames,11)

msp = SpeechBox.specgram(frames)
pp = SpeechBox.pitch_spec_comb(frames)

