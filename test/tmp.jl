using WAV
using Plots
using BenchmarkTools
# using DSP

testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"))
x, fs = wavread(joinpath(testdir,"King.wav"))

frames = framed_signal(x,fs,0.09,0.01)
frame = extract_frame(frames,11)

h = SpeechBox.lpc_response(frame,fs)

lpspec = SpeechBox.lpc_response(frames)

mag = SpeechBox.magspec(frame,fs)

@btime pe = SpeechBox.generate_pitch_estimate(frame,fs)

# msp = specgram(frames)

# psp = periodogram(frames)

# plot(log(msp))
