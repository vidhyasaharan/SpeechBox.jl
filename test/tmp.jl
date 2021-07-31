using SpeechBox
using WAV
using Plots; plotlyjs()
using BenchmarkTools
# using TIMITutilities

testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"))
x, fs = wavread(joinpath(testdir,"SA1.wav"))

# utt = TIMITutilities.uttlist(TRAIN)
# x, fs = TIMITutilities.readriff(utt[2])

s = speech_waveform(x,fs)

frames = framed_signal(s,0.05,0.01)
frame = extract_frame(frames,32)

dd = zeros(Float64,frames.num_signal_frames-1)
win = SpeechBox.window(frames.frame_length; wtype="hamming")
p = SpeechBox.lpc_order(frames.signal.fs)

for i ∈ eachindex(dd)
    f1 = extract_frame(frames,i)
    f2 = extract_frame(frames,i+1)
    dd[i] = SpeechBox.distitak2(f1.*win,f2.*win,p)
end




ptf = SpeechBox.RAPT_timefreq(s)
ca = SpeechBox.RAPT_pitch_candidates(s)
lc = SpeechBox.RAPT_cands_local_costs(ca)


# msp = SpeechBox.specgram(frames)
# pp = SpeechBox.pitch_spec_comb(frames)

# cr, frqs = SpeechBox.xcorr_spectral_comb(frames)
# h,z = SpeechBox.generate_logfrq_pitch_comb(;frq_per_octave=200)
# pd = SpeechBox.periodogram(frames,frqs)
# lpd = log.(pd.components)

# fnum = 145
# ii,mm = SpeechBox.findpeaks_sorted(cr[:,fnum])
# p = plot(cr[:,fnum])
# plot!(p,ii,mm, seriestype=:scatter)