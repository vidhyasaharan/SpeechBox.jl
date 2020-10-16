using SpeechBox
using WAV

testdir = normpath(joinpath(dirname(@__FILE__),"../test/"))
srcdir = normpath(joinpath(dirname(@__FILE__),"../src/"))

x, fs = wavread(joinpath(testdir,"King.wav"))
x = x[:]
signal = speech_waveform(x,fs)
sig_frames = framed_signal(signal)
sframe = extract_frame(sig_frames,52)

frqs = SpeechBox.logfreq_array(fmin = 1, fmax = fs/2, frq_per_octave = 1000)
pd = periodogram(sframe,fs,frqs)
