using WAV
using BenchmarkTools
using Plots

testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"))
x, fs = wavread(joinpath(testdir,"King.wav"))

frames = framed_signal(x,fs,0.09,0.01)
frame = extract_frame(frames,11)

ce = SpeechBox.cexp(100,fs,length(frame))

@benchmark t1 = abs2(SpeechBox.dot_noavx(frame,ce))
@benchmark t2 = abs2(SpeechBox.dotavx(frame,ce))


fmin = 10
fmax = fs/2
frqs = SpeechBox.logfreq_array(;fmin = fmin,fmax = fmax)


@benchmark pd = SpeechBox.periodogram(frame,fs,frqs)
@benchmark pd1 = SpeechBox.periodogram_avx(frame,fs,frqs)
@benchmark pd2 = SpeechBox.periodogram_mul(frame,fs,frqs)

@benchmark pd = SpeechBox.periodogram(frames,frqs)
@benchmark pd1 = SpeechBox.periodogram_avx(frames,frqs)