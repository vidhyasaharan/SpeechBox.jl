using Revise
using SpeechBox
using Makie
using WAV


testdir = normpath(joinpath(dirname(@__FILE__),"../test/"))
srcdir = normpath(joinpath(dirname(@__FILE__),"../src/"))

x, fs = wavread(joinpath(testdir,"King.wav"))


function plot_freqz(filter,fs::Number)
    fz = DSP.freqz(filter,0:1:fs/2,fs)
    sc = lines(log.(abs.(fz)))
    return sc
end


n_res = 5;
fs = 8000;

ff,frqs,bws = SpeechBox.rand_allpole(fs,n_res)

iz = DSP.impz(ff,500)
nz = DSP.filt(ff,randn(500))


f1 = SpeechBox.frame_lpc(iz,10)
f2 = SpeechBox.frame_lpc(nz,10)

fz = DSP.freqz(ff,0:1:fs/2,fs)
fz1 = DSP.freqz(f1,0:1:fs/2,fs)
fz2 = DSP.freqz(f2,0:1:fs/2,fs)

sc = Scene()
lines!(sc,log.(abs.(fz)))
lines!(sc,log.(abs.(fz2)),color=:red)





# win_dur = 0.03
# win_overlap = 0.01
# sig_frames = framed_signal(x,win_dur,win_overlap)
# vi1 = vad(sig_frames,alg = "energy_threshold",energy_threshold = 0.1)
# vi2 = vad(sig_frames,alg = "unvoiced_fraction",unvoiced_fraction = 0.4)
#
# vi = vad_energy_threshold(sig_frames,0.05)
# vi2 = vad_energy_fraction(sig_frames,0.3)
#
# en = frame_energy(sig_frames,normalised=true)
#
# scene = lines(vi,color=:blue)
# lines!(scene,en,color=:red)
# lines!(scene,vi1,color=:green)
# lines!(scene,vi2,color=:black)
