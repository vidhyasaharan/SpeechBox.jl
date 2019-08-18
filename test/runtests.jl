
using SpeechBox
using Test
using LibSndFile
using FileIO
using SampledSignals


testdir = normpath(joinpath(dirname(@__FILE__),"../test/"))
srcdir = normpath(joinpath(dirname(@__FILE__),"../src/"))

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

@testset "magspec" begin
    t = 0:0.01:0.99
    frq = 25
    xx = cos.(2*pi*frq*t)
    mspec = magspec(xx,100.0;wtype = "rect")
    mmag,mfrq = findmax(mspec)
    @test length(mspec) == 51
    @test mfrq == frq + 1
    @test round(mmag) == 50.0
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

@testset "melfcc" begin
    ncof = 14
    nfil = 22
    win_dr = 0.025
    win_ovlp = 0.01
    mfc = melfcc(x;ncoef = ncof, nfilt = nfil, win_dur = win_dr, win_overlap = win_ovlp)
    frames = framed_signal(x,win_dr,win_ovlp)
    fx = extract_frame(frames,5)
    @test size(mfc,1) == ncof
    @test size(mfc,2) == frames.num_frames
end
