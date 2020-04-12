#Generate window function of given length
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

#Generate array of frequencies (in Hz), equally spaced in log domain with resolution given in frequencies per octave
function logfreq_array(;fmin::Number = 10, fmax::Number = 4000, frq_per_octave::Number = 120)
    lfmin = log2(fmin)
    lfmax = log2(fmax)
    lfres = 1/frq_per_octave
    lfrq = lfmin:lfres:lfmax
    return exp2.(lfrq)
end

#Generate complex exponential sequence
cexp(f,fs,N) = (1/sqrt(N))*exp.(2π*im*(f/fs)*(1:N))

#Generate projection matrix for complex exponential signals/vectors
function periodogram_basis_matrix(frqs::Array{T,1},fs::Number,N::Int) where T<:Number
    nfrqs = length(frqs)
    periodogram_basis = zeros(Complex{Float},nfrqs,N)
    for i=1:nfrqs
        periodogram_basis[i,:] = cexp(frqs[i],fs,N)
    end
    return periodogram_basis
end
