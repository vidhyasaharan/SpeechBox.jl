
function vad_energy_threshold(sig_frames::framed_signal, energy_thr::Float = 0.05)
    energy = frame_energy(sig_frames)
    max_energy = maximum(energy)
    abs_thr = max_energy*energy_thr
    v_indx = zeros(Int,sig_frames.num_frames)
    for i in findall(energy.>abs_thr)
        setindex!(v_indx,1,i)
    end
    return v_indx
end

function vad_energy_fraction(sig_frames::framed_signal, unvoiced_fraction::Float = 0.1)
    energy = frame_energy(sig_frames)

    sorted_energy = sort(energy)
    num_unvoiced_frames = Int(round(unvoiced_fraction*sig_frames.num_frames))
    en_thr = sorted_energy[num_unvoiced_frames]

    v_indx = zeros(Int,sig_frames.num_frames)
    for i in findall(energy.>en_thr)
        setindex!(v_indx,1,i)
    end
    return v_indx
end

function vad(sig_frames::framed_signal;alg = "energy_threshold", params...)
    if alg == "energy_threshold"
        if haskey(params,:energy_threshold)
            en_thr = params[:energy_threshold]
        else
            println("Warning: energy_threshold parameter value not supplied, using default value of 0.05")
            en_thr = 0.05
        end
        return vad_energy_threshold(sig_frames,en_thr)
    elseif alg == "unvoiced_fraction"
        if haskey(params, :unvoiced_fraction)
            uv_fr = params[:unvoiced_fraction]
        else
            println("Warning: unvoiced_fraction parameter value not supplied, using default value of 0.2")
            uv_fr = 0.2
        end
        return vad_energy_fraction(sig_frames,uv_fr)
    else
        println("Error: Unrecognised VAD algorithm")
        return 0
    end
end
