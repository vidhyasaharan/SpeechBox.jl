using Revise

testdir = normpath(joinpath(dirname(@__FILE__),"../test/"))
srcdir = normpath(joinpath(dirname(@__FILE__),"../src/"))

# includet(joinpath(srcdir,"Speechbox.jl"))
# using .Speechbox
using SpeechBox
using Test
using LibSndFile
using FileIO
using Makie
using SampledSignals

# using Plots

x = load(joinpath(testdir,"King.wav"))
fs = samplerate(x)



win_dur = 0.03
win_overlap = 0.01
sig_frames = framed_signal(x,win_dur,win_overlap)
vi1 = vad(sig_frames,alg = "energy_threshold",energy_threshold = 0.1)
vi2 = vad(sig_frames,alg = "unvoiced_fraction",unvoiced_fraction = 0.4)

vi = vad_energy_threshold(sig_frames,0.05)
vi2 = vad_energy_fraction(sig_frames,0.3)

en = frame_energy(sig_frames,normalised=true)

scene = lines(vi,color=:blue)
lines!(scene,en,color=:red)
lines!(scene,vi1,color=:green)
lines!(scene,vi2,color=:black)
