
#Convert analogue frequency (f) in Hz to digital frequency (θ) given fs
freq2θ(f::Real, fs::Real) = 2pi*(f/fs)

function freq2θ(f::AbstractVector{<:Real}, fs::Real)
    θ = Vector{Float}(undef,length(f))
    freq2θ!(θ,f,fs)
    return θ
end

function freq2θ!(θ::AbstractVector{<:Real}, f::AbstractVector{<:Real}, fs::Real)
    k = 2pi/fs
    for i ∈ eachindex(θ)
        θ[i] = k*f[i]
    end
    return θ
end

#Map analogue frequency (f) in Hz to unit circle on z-plane given fs
freq2z(f::Real,fs::Real) = exp(im*freq2θ(f,fs))

function freq2z(f::AbstractVector{<:Real}, fs::Real)
    z = Vector{Complex{Float}}(undef,length(f))
    freq2z!(z,f,fs)
    return z
end

function freq2z!(z::AbstractVector{<:Complex}, f::AbstractVector{<:Real}, fs::Real)
    for i ∈ eachindex(z)
        z[i] = exp(im*freq2θ(f[i],fs))
    end
end

#Filter structs
struct filter_coefs
    num::Array{Float}
    den::Array{Float}
end


function powers(x::Number, N::Int)
    px = Vector{eltype(x)}(undef,N+1)
    px[1] = one(eltype(px))
    for i ∈ 1:N
        px[i+1] = x^i
    end
    return px
end


function H(F::filter_coefs, z::Complex{<:Real})
    npz = powers(z, length(F.num)-1)
    dpz = powers(z, length(F.den)-1)
    Hz = dot(F.num,npz)/dot(F.den,dpz)
    return Hz
end

H(F::filter_coefs, f::Real, fs::Real) = H(F,freq2z(f,fs))


function filter_resp(F::filter_coefs, f::AbstractVector{<:Real}, fs::Real)
    z = freq2z(f,fs)
    Hz = map(x->H(F,x), z)
    return Hz
end


function Hmag(F::filter_coefs, θ::Real)
    iθ = im*θ
    Nm = convert(Complex{Float},F.num[1])
    Dm = convert(Complex{Float},F.den[1])
    Nord = length(F.num) - 1
    Dord = length(F.den) - 1
    if(Nord>0)
        for i ∈ 1:Nord
            Nm += F.num[i+1]*exp(i*iθ)
        end
    end
    if(Dord>0)
        for i ∈ 1:Dord
            Dm += F.den[i+1]*exp(i*iθ)
        end
    end
    return abs(Nm)/abs(Dm)
end

Hmag(F::filter_coefs, f::Real, fs::Real) = Hmag(F, freq2θ(f,fs))


filter_magresp(F::filter_coefs, f::AbstractVector{<:Real}, fs::Real) = map(x->Hmag(F,x,fs),f)