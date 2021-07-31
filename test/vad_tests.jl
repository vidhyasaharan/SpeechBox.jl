@testset "vad_energy" begin
    win_dur = 0.03
    win_shift = 0.01
    sig_frames = framed_signal(x,fs,win_dur,win_shift)

    en_thr = 0.05
    # vi = vad(sig_frames,alg = "energy_threshold",energy_threshold = en_thr)
    vi = vad(energy_threshold(), sig_frames,threshold = en_thr)


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
