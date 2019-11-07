
using Test
using SpeechBox
# using LibSndFile
# using FileIO
# using SampledSignals
using WAV

testdir = normpath(joinpath(dirname(@__FILE__),"../test/"))
srcdir = normpath(joinpath(dirname(@__FILE__),"../src/"))

x, fs = wavread(joinpath(testdir,"King.wav"))

# x = load(joinpath(testdir,"King.wav"))
# fs = samplerate(x)

@testset "framed_signal" begin
    # frames = framed_signal(x) #For SampleBuf input
    # @test typeof(frames) == framed_signal
    # @test typeof(frames.x) <: Array{Float}
    # @test nframes(frames.x) > 0
    # @test typeof(frames.frame_length) <: Int
    # @test typeof(frames.frame_overlap) <: Int
    # @test typeof(frames.num_signal_frames) <: Int
    # @test typeof(frames.num_frames) <: Int
    # @test frames.frame_length > 0
    # @test frames.frame_overlap > 0
    # @test frames.num_signal_frames > 0
    # @test frames.num_frames > 0
    # @test frames.num_frames >= frames.num_signal_frames
    # @test frames.num_signal_frames >= (nframes(frames.x)-frames.frame_length)/frames.frame_overlap

    frames = framed_signal(x,fs) #for Array{Float} input
    @test typeof(frames) == framed_signal
    @test typeof(frames.x) <: Array{Float}
    @test length(frames.x) > 0
    @test typeof(frames.frame_length) <: Int
    @test typeof(frames.frame_overlap) <: Int
    @test typeof(frames.num_signal_frames) <: Int
    @test typeof(frames.num_frames) <: Int
    @test frames.frame_length > 0
    @test frames.frame_overlap > 0
    @test frames.num_signal_frames > 0
    @test frames.num_frames > 0
    @test frames.num_frames >= frames.num_signal_frames
    @test frames.num_signal_frames >= (length(frames.x)-frames.frame_length)/frames.frame_overlap

end

@testset "extract_frame" begin
    frames = framed_signal(x,fs)
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
    msp = specgram(x, fs;win_dur = 0.03,win_overlap=0.01, wtype = "hamming")
    @test typeof(msp) <: Array{<:AbstractFloat,2}
    @test size(msp,1) > 0
    @test size(msp,2) > 0
    @test maximum(isa.(msp,Complex))==false
    @test minimum(msp) >= eps()

    t = 0:0.01:0.99
    frq1 = 25
    frq2 = 30
    xx = [cos.(2*pi*frq1*t);cos.(2*pi*frq2*t)]
    msp = specgram(xx,100.0,win_dur = 1.0, win_overlap = 1.0)
    @test size(msp,1) > 0
    @test size(msp,2) > 0
    @test maximum(isa.(msp,Complex))==false
    @test minimum(msp) >= eps()
    mmag1,mfrq1 = findmax(msp[:,1])
    mmag2,mfrq2 = findmax(msp[:,2])
    @test mfrq1 == frq1 + 1
    @test mfrq2 == frq2 + 1

end

@testset "melbankm" begin
    mb = melbankm(fs,Int(ceil(0.02*fs)),nfilt = 22)
    @test size(mb,1) == 22
    @test size(mb,2) == Int(ceil(0.02*fs))
    @test sum(mb,dims=1) == ones(1,Int(ceil(0.02*fs)))
end

@testset "vad_energy" begin
    win_dur = 0.03
    win_overlap = 0.01
    sig_frames = framed_signal(x,fs,win_dur,win_overlap)

    en_thr = 0.05
    vi = vad(sig_frames,alg = "energy_threshold",energy_threshold = en_thr)


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

@testset "melfcc" begin
    ncof = 14
    nfil = 22
    win_dr = 0.025
    win_ovlp = 0.01
    frames = framed_signal(x,fs,win_dr,win_ovlp)
    mfc = melfcc(x,fs;ncoef = ncof, nfilt = nfil, win_dur = win_dr, win_overlap = win_ovlp)
    fx = extract_frame(frames,5)
    @test size(mfc,1) == ncof
    @test size(mfc,2) == frames.num_frames
end
