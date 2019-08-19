function window(flen::Int;wtype::String="hanning")
    if(wtype=="rect")
        win = ones(flen)
    elseif(wtype=="hamming")
        win = hamming(flen)
    elseif(wtype=="hanning")
        win = hanning(flen)
    else
        println("Warning: window type not recognised - using Hann window")
        win = hanning(flen)
    end
    return win
end
