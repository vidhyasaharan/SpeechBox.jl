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


@testset "framed_signal" begin
    frames = framed_signal(x)
    @test typeof(frames) == framed_signal
    @test typeof(frames.x) <: AbstractSampleBuf
    @test nframes(frames.x) > 0
    @test typeof(frames.frame_length) <: Int
    @test typeof(frames.frame_overlap) <: Int
    @test typeof(frames.num_signal_frames) <: Int
    @test typeof(frames.num_frames) <: Int
    @test frames.frame_length > 0
    @test frames.frame_overlap > 0
    @test frames.num_signal_frames > 0
    @test frames.num_frames > 0
    @test frames.num_frames >= frames.num_signal_frames
    @test frames.num_signal_frames >= (nframes(frames.x)-frames.frame_length)/frames.frame_overlap
end

@testset "extract_frame" begin
    frames = framed_signal(x)
    @test typeof(extract_frame(frames,1))<:Array{<:AbstractFloat,1}
    @test length(extract_frame(frames,1)) > 0
    @test length(extract_frame(frames,frames.num_frames)) > 0
    @test maximum(abs.(extract_frame(frames,frames.num_frames))) > 0
end

@testset "specgram" begin
    msp = specgram(x;win_dur = 0.03,win_overlap=0.01, wtype = "hamming")
    @test typeof(msp) <: Array{<:AbstractFloat,2}
    @test size(msp,1) > 0
    @test size(msp,2) > 0
    @test maximum(isa.(msp,Complex))==false
    @test minimum(msp) >= eps()
end

@testset "melbankm" begin
    mb = melbankm(fs,Int(ceil(0.02*fs)),nfilt = 22)
    @test size(mb,1) == 22
    @test size(mb,2) == Int(ceil(0.02*fs))
    @test sum(mb,dims=1) == ones(1,Int(ceil(0.02*fs)))
end

@testset "vad_energy" begin
    using DSP

    win_dur = 0.03
    win_overlap = 0.01
    sig_frames = framed_signal(x,win_dur,win_overlap)

    en_thr = 0.05
    vi = vad_energy_threshold(sig_frames,en_thr)

    energy = zeros(Float,sig_frames.num_frames)
    for i=1:sig_frames.num_frames
        frame = extract_frame(sig_frames,i)
        energy[i] = rms(frame);
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


# msp1 = specgram(x)
# msp2 = specgram(x,wtype="rect")
# msp3 = specgram(x,wtype="hamming")
#
# image(msp1')




vi = vad_energy_threshold(sig_frames,0.05)

scene = lines(vi,color=:blue)
lines!(scene,en,color=:red)
