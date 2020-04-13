using WAV

testdir = normpath(joinpath(dirname(@__FILE__),"../test/"))
x, fs = wavread(joinpath(testdir,"King.wav"))

frames = framed_signal(x,fs,0.09,0.01)
frame = extract_frame(frames,11)

fsig = speech_waveform(frame,fs)

t = (0:length(fsig.x)-1)/fs


msp = specgram(frames)
nfrqs = size(msp,1)
frqs = convert.(Float,collect(range(0, fs/2, length = nfrqs)))

mspec = magspec(frame,fs)
pgram = periodogram(frame,fs)
