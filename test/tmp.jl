using WAV
using Plots; plotly()

testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"))
x, fs = wavread(joinpath(testdir,"piano2.wav"))

frames = framed_signal(x,fs,0.09,0.01)
frame = extract_frame(frames,11)

mag = SpeechBox.magspec(frame,fs)

msp = specgram(frames)

# psp = periodogram(frames)

plot(log(msp))
