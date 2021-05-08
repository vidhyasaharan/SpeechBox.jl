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
v = SpeechBox.vad(frames)
pin = SpeechBox.frqindex(pp,msp.frqs)

p = plot(log(msp))
plot!(p,1:frames.num_signal_frames,pin, c = :green, legend = false)
plot!(p,[15],seriestype = :vline, c = :green, legend = false)
plot!(log(msp))
plot!(p,1:frames.num_signal_frames,200*ones(frames.num_signal_frames),c = :green)
# psp = periodogram(frames)

# plot(log(msp))
