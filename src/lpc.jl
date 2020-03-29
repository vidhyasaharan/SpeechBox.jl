function allpole(pf::Array{<:Number,1} = [1000, 1800, 2900, 3400, 5000, 6800], pbw::Array{<:Number,1} = [50, 120, 200, 300, 500, 800], fs::Number = 16000)
    num_poles = length(pf)
    poles  = zeros(Complex{Float},2*num_poles)
    for i=1:num_poles
        r = exp(-pi*pbw[i]/fs)
        # println(r)
        poles[2i-1] = r*exp(1im*2*pi*pf[i]/fs)
        poles[2i] = conj(poles[2i-1])
    end
    return DSP.Filters.ZeroPoleGain([0],poles,1)
end

function rand_vocalfilter_resonance(fmin::Number = 0, fmax::Number = 8000)
    frequency = fmin + ((fmax-fmin)*rand())
    bw_min = 0.08*frequency
    bw_max = 0.2*frequency
    bandwidth = bw_min + ((bw_max-bw_min)*rand())
    return frequency, bandwidth
end

function rand_allpole(fs::Number = 8000, num_res::Number = 10)
    fmax = fs/2
    fint = fmax/num_res
    frqs = zeros(Float,num_res)
    bws = zeros(Float,num_res)
    flo = 50
    fhi = 1.5*fint
    for i=1:num_res
        frqs[i],bws[i] = rand_vocalfilter_resonance(flo, fhi)
        fmid = i*fint + 0.5*fint
        flo = fmid-fint
        if(frqs[i]>flo)
            flo = frqs[i]
        end
        fhi = fmid+fint
    end
    ar_filt = allpole(frqs,bws,fs)
    return ar_filt, frqs, bws
end


function frame_lpc(x::Array{Float,1},order::Int)
    a = ones(Float,order+1)
    temp,err = DSP.LPC.lpc(x,order)
    a[2:end] = temp
    filter = DSP.Filters.PolynomialRatio([1],a)
    return filter
end
