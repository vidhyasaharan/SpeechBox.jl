using WAV

testdir = normpath(joinpath(dirname(@__FILE__),"../test/"))
x, fs = wavread(joinpath(testdir,"King.wav"))

frames = framed_signal(x,fs,0.09,0.01)
frame = extract_frame(frames,11)


msp = specgram(frames)
