"""
    vad_energy_threshold(sig_frames::framed_signal, energy_thr=0.05)

Estimates the signal 'energy' in each frame (L2 norm computed using the `frame_energy` function) and assigns every frame with 'energy' greater than `energy_thr` times the maximum frame energy as voiced.
"""
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

"""
    vad_energy_fraction(sig_frames::framed_signal, unvoiced_fraction=0.2)

Estimates the signal 'energy' in each frame (L2 norm computed using the `frame_energy` function) and assigns a certain fraction of the frames (specified by `unvoiced_fractions`) with the lowest energy as not voiced.
"""
function vad_energy_fraction(sig_frames::framed_signal, unvoiced_fraction::Float = 0.2)
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


"""
    vad(sig_frames::framed_signal; <keyword arguments>)

Determine if each frame in the framed\\_signal object `sig_frames` corresponds to voiced speech or not
# Arguments
- `alg` : the vad algorithm to use. The options are: (i) `energy_threshold` [default] which compares normalised frame energy to a fixed threshold (keyword argment); or (ii) `unvoiced_fraction` which sets a certain fraction (keyword argument) of lowest energy frames as not voiced
- `energy_threshold` : minimum frame energy for a frame to be considered voiced when using `alg = "energy_threshold"`  [Default = 0.05]
- `unvoiced_fraction` : fraction of frames that will be marked as not voiced with using `slg = "unvoiced_fraction"` [Default = 0.2]
"""
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
