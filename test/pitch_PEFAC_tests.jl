using SpeechBox
using WAV

testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"))
srcdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../src/"))

x, fs = wavread(joinpath(testdir,"King.wav"))
x = x[:]
signal = speech_waveform(x,fs)
sig_frames = framed_signal(signal)
sframe = extract_frame(sig_frames,52)


frqs = SpeechBox.logfreq_array(fmin = 1, fmax = fs/2, frq_per_octave = 1000)

ltass_mag,frqs = SpeechBox.LTASS(fs/2,1000)
pd = SpeechBox.periodogram(sig_frames,frqs)

lpd = pow2db(pd)



npd = SpeechBox.LTASS_normalise(pd)
plot(npd)
plot(npd.components[:,149])


p1 = npd.components[:,52]
frq_per_octave = 1000
γ = 1.8
K = 5
h,z = SpeechBox.generate_logfrq_pitch_comb(γ,K, frq_per_octave = frq_per_octave)
padded_pd = zeros(length(p1)+length(h)-1)
padded_pd[z:z+length(p1)-1] = p1
y = zeros(size(p1))
for i=1:length(p1)
    y[i] = dot(padded_pd[i:i+length(h)-1],h)
end

tt = SpeechBox.generate_pitch_estimate(sig_frames)
