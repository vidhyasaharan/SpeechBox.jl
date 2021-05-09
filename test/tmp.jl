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

cr, frqs = SpeechBox.xcorr_spectral_comb(frames)
h,z = SpeechBox.generate_logfrq_pitch_comb(;frq_per_octave=200)
pd = SpeechBox.periodogram(frames,frqs)
lpd = log.(pd.components)

fnum = 145
ii,mm = SpeechBox.findpeaks_sorted(cr[:,fnum])
p = plot(cr[:,fnum])
plot!(p,ii,mm, seriestype=:scatter)