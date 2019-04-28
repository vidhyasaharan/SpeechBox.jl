using Revise

testdir = normpath(joinpath(dirname(@__FILE__),"../test/"))
srcdir = normpath(joinpath(dirname(@__FILE__),"../src/"))

# includet(joinpath(srcdir,"Speechbox.jl"))
# using .Speechbox
using SpeechBox
using Test
using LibSndFile
using FileIO
using Makie
using SampledSignals

# using Plots

x = load(joinpath(testdir,"King.wav"))
fs = samplerate(x)


@testset "vad_energy" begin
    using DSP

    win_dur = 0.03
    win_overlap = 0.01
    sig_frames = framed_signal(x,win_dur,win_overlap)

    en_thr = 0.05
    vi  = vad_energy_threshold(sig_frames,en_thr)

    energy = zeros(Float,sig_frames.num_frames)
    for i=1:sig_frames.num_frames
        frame = extract_frame(sig_frames,i)
        energy[i] = sum(abs2,frame)
        # energy[i] = rms(frame);
    end
    max_energy = maximum(energy)
    energy = energy./max_energy

    for i=1:sig_frames.num_frames
        if vi[i] == 1
            @test energy[i] > en_thr
        else
            @test energy[i] <= en_thr
        end
    end
end

win_dur = 0.03
win_overlap = 0.01
sig_frames = framed_signal(x,win_dur,win_overlap)
vi = vad_energy_threshold(sig_frames,0.05)
vi2 = vad_energy_fraction(sig_frames,0.3)

en = frame_energy(sig_frames)

scene = lines(vi,color=:blue)
lines!(scene,en/50,color=:red)
lines!(scene,vi2,color=:green)
