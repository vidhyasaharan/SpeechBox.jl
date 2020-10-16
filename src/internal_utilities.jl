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
cexp(f,fs,N) = exp.(2π*im*(f/fs)*(1:N))

#Generate projection matrix for complex exponential signals/vectors
function cexp_proj_matrix(frqs::Array{T,1},fs::Number,N::Int) where T<:Number
    nfrqs = length(frqs)
    proj_matrix = zeros(Complex{Float},nfrqs,N)
    for i=1:nfrqs
        proj_matrix[i,:] = cexp(frqs[i],fs,N)
    end
    return proj_matrix
end


#Generate first order difference of a sequence y[i] = x[i+1] - x[i] (output sequence length is 1 less than input sequence length)
function Δ(x::Vector)
    len = length(x)
    Δx = zeros(typeof(x[1]),len-1)
    if(len>1)
        for i=1:len-1
            Δx[i] = x[i+1]-x[i]
        end
    end
    return Δx
end


#Find local peaks/maximas in a sequence
function findpeaks(x::Vector; min_dist::Int = 2)
    ind = Int[]
    mag = Real[]
    if x[1]>x[2]
        push!(ind,1)
        push!(mag,x[1])
    end
    for i=2:length(x)-1
        if(x[i-1]<x[i]>x[i+1])
            push!(ind,i)
            push!(mag,x[i])
        end
    end
    if(x[end]>x[end-1])
        push!(ind,length(x))
        push!(mag,x[end])
    end
    while(minimum(Δ(ind))<min_dist)
        ind,mag = remove_nearest_peak(ind,mag)
    end
    return ind,mag
end


#Support function for findpeaks() - removes the smaller of the two closest peaks in a set of local peaks
function remove_nearest_peak(ind::Vector,mag::Vector)
    npks = length(ind)
    if(npks>1)
        dist = Δ(ind)
        m_i = argmin(dist)
        if(mag[m_i+1]<mag[m_i])
            m_i += 1
        end
        deleteat!(ind,m_i)
        deleteat!(mag,m_i)
    end
    return ind,mag
end


function resample(signal::speech_waveform, fs_new::Number)
    rx = DSP.Filters.resample(signal.x, fs_new/signal.fs)
    fs = convert(Float,fs_new)
    return speech_waveform(rx,fs)
end



function element_op_spectrum(func::AbstractString)
    me = Expr(:call, :map, Meta.parse(func), :(sp.components))
    re = :(spectrum(sp.signal,$me,sp.frqs,sp.title))
    le = Expr(:call, Meta.parse(func), :(sp::SpeechBox.spectrum))
    return Expr(:(=), le, re)
end

function element_op_timefreq(func::AbstractString)
    me = Expr(:call, :map, Meta.parse(func), :(tf.components))
    re = :(timefreq(tf.signal,tf.frames,$me,tf.frqs,tf.time,tf.title))
    le = Expr(:call, Meta.parse(func), :(tf::SpeechBox.timefreq))
    return Expr(:(=), le, re)
end

eval(element_op_spectrum("Base.log10"))
eval(element_op_spectrum("Base.log"))
eval(element_op_spectrum("DSP.amp2db"))
eval(element_op_spectrum("DSP.pow2db"))

eval(element_op_timefreq("Base.log10"))
eval(element_op_timefreq("Base.log"))
eval(element_op_timefreq("DSP.amp2db"))
eval(element_op_timefreq("DSP.pow2db"))
