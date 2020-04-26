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

# function Base.log(tf::timefreq)
#     return timefreq(tf.signal,tf.frames,log.(tf.components),tf.frqs,tf.time,tf.title)
# end
#
# function Base.log(sp::spectrum)
#     return spectrum(sp.signal,log.(sp.components),sp.frqs,sp.title)
# end

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

eval(element_op_timefreq("Base.log10"))
eval(element_op_timefreq("Base.log"))
