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
    fmin = Float(fmin)::Float
    fmax = Float(fmax)::Float
    frq_per_octave = Float(frq_per_octave)::Float
    lfmin = log2(fmin)
    lfmax = log2(fmax)
    lfres = 1/frq_per_octave
    lfrq = lfmin:lfres:lfmax
    return exp2.(lfrq)
end

#Generate array of desired number of equally spaced frequencies (in Hz)
function linfreq_array(;fmin::Number = 0, fmax::Number = 4000, nfrqs::Int = 80)
    fmin = Float(fmin)::Float
    fmax = Float(fmax)::Float
    fres = (fmax-fmin)/(nfrqs-1)
    frqs = fmin:fres:fmax
    return collect(frqs)
end



#Generate complex negative exponential sequence with norm = 1
cexp(f::T,fs::T,N::Int) where {T<:Real} = (1/sqrt(N))*exp.(-2π*im*(f/fs)*(1:N))


#Generate projection matrix for complex exponential signals/vectors
function cexp_proj_matrix(frqs::Vector{<:Real},fs::Real,N::Int)
    nfrqs = length(frqs)
    proj_matrix = zeros(Complex{Float},nfrqs,N)
    for i in eachindex(frqs)
        proj_matrix[i,:] = cexp(frqs[i],fs,N)
    end
    return proj_matrix
end


#Generate first order difference of a sequence y[i] = x[i+1] - x[i] (output sequence length is 1 less than input sequence length)
function Δ(x::AbstractVector)
    len = length(x)
    Δx = zeros(typeof(x[1]),len-1)
    if(len>1)
        @views for i=1:len-1
            Δx[i] = x[i+1]-x[i]
        end
    end
    return Δx
end


#Find local peaks/maximas in a sequence
function findpeaks(x::Vector; min_dist::Int = 2)
    ind = Int[]
    mag = eltype(x)[]
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
    if(length(ind)>1)
        while(minimum(Δ(ind))<min_dist)
            remove_nearest_peak!(ind,mag)
        end
    end
    return ind,mag
end



#Local peaks/maximas sorted by magnitude
function findpeaks_sorted(x::Vector; min_dist::Int = 2, num_peaks::Int = 0)
    ind, mag = findpeaks(x;min_dist)
    si = sortperm(mag, rev=true)
    if((num_peaks>0)&&(num_peaks<length(mag)))
        si = si[1:num_peaks]
    end
    return ind[si], mag[si]
end



#Support function for findpeaks() - removes the smaller of the two closest peaks in a set of local peaks
function remove_nearest_peak!(ind::Vector,mag::Vector)
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
    # return ind,mag
end

#Resample signal in speech_waveform object (wrapper for resample from DSP.jl)
function resample(signal::speech_waveform, fs_new::Number)
    rx = DSP.Filters.resample(signal.x, fs_new/signal.fs)
    fs = convert(Float,fs_new)
    return speech_waveform(rx,fs)
end

#Index of closest frequency in an array to a given frequency
frqindex(f::Float, frqs::AbstractVector{Float}) = argmin(abs.(frqs.-f))

function frqindex(f::AbstractVector{Float}, frqs::AbstractVector{Float})
    findx = Vector{Int}(undef,length(f))
    for i in eachindex(f)
        findx[i] = frqindex(f[i],frqs)
    end
    return findx
end

#Sample index closest in time to a given time (in secs) based on sampling frequency
timeindex(t::Real, fs::Real) = Int(round(t*fs))


#Zero pad a vector
zero_pad(x::Vector, pad_len::Int) = [zeros(pad_len); x; zeros(pad_len)]

#Symmetrically pad a vector
symmetric_pad(x::Vector, pad_len::Int) = [x[pad_len+1:-1:2]; x; x[end-1:-1:end-(pad_len)]]

#Remove padding from a vector
unpad_vector(x::Vector, pad_len::Int) = x[pad_len+1:end-pad_len]

#Generate white noise
white_noise(len::Int) = randn(MersenneTwister(), Float, len)
white_noise(dur::Number, fs::Number) = white_noise(dur2len(dur,fs))

#Generate AR process noise
ar_process(a::Vector, len::Int) = DSP.filt([1],a, white_noise(len))
ar_process(a::Vector, dur::Number, fs::Number) = ar_process(a,dur2len(dur,fs))

#Infer number of samples from duration and sampling rate
dur2len(dur::Number, fs::Number) = Int(round(dur*fs))

#Generate impulse train
function impulse_train(period::Int, len::Int)
    x = zeros(Float,len)
    x[1:period:end] .= 1
    return x
end

impulse_train(f₀::Real, dur::Real, fs::Real) = impulse_train(Int(round(fs/f₀)), Int(round(dur*fs)))


#Elementwise operations on spectrum and timefreq components
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
