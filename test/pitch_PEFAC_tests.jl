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

ltass_mag,frqs = SpeechBox.LTASS(fs/2,1000)
pd = SpeechBox.periodogram(sig_frames,frqs)
sm = SpeechBox.moving_average(pow2db(pd),time_window = 11, freq_window = 11)
lpd = pow2db(pd)


fnum = 55
tt = lpd.components[:,fnum].*ltass_mag./sm.components[:,fnum]
plot(tt)

plot(sm)
plot(lpd)


npd = SpeechBox.LTASS_normalise(lpd)
plot(npd)
plot(npd.components[:,149])

tt = maximum(lpd.components)
