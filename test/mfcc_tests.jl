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
    frames = framed_signal(x,fs,win_dr,win_ovlp)
    mfc = melfcc(x,fs;ncoef = ncof, nfilt = nfil, win_dur = win_dr, win_shift = win_ovlp)
    fx = extract_frame(frames,5)
    @test size(mfc,1) == ncof
    @test size(mfc,2) == frames.num_frames
end
