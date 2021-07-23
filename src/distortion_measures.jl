#Itakura Saito Distortion
function distispf(x::AbstractVector{Float}, y::AbstractVector{Float}, p::Int)
    npts = 1000
    Δθ = 2pi/(npts-1)
    θ = -pi:Δθ:pi
    pf1 = lpc_magz(x,2pi,p; frqs = θ)
    pf2 = lpc_magz(y,2pi,p; frqs = θ)
    d = zero(Float)
    for i ∈ eachindex(pf1)
        v = pf1[i]/pf2[i]
        t = v - log(v) - 1
        d += t
    end
    d *= Δθ/2pi
    return d
end


function distis(x::AbstractVector{Float}, y::AbstractVector{Float}, p::Int)
    a = lpc(x, p)
    b = lpcar2ra(a)
    for i ∈ 2:length(b)
        b[i] *= 2
    end

    Ryy = acorr(y, p+1)
    â = lpc(Ryy)
    r = lpcacorr2v(Ryy)

    c = log(dotavx(a))
    # d = c + log(dotavx(b,r)) - log(dotavx(â,r))
    d = dotavx(b,r)/dotavx(â,r)
    return d
end


function distis_mat(x::AbstractVector{Float}, y::AbstractVector{Float}, p::Int)
    rxx = acorr(x, p+1)
    V = acorr_mat(rxx)
    a = lpc(x, p)
    â = lpc(y, p)
    N = dotavx(a, V*a)
    D = dotavx(â, V*â)
    d = N/D
    return d
end


#Convert inverse filter coefficients to autocorrelation coefficients
function lpcar2ra(a::AbstractVector{Float})
    na = 1/dotavx(a)
    p = length(a)
    b = Vector{Float}(undef,p)
    b[1] = 1
    for i ∈ 2:p
        b[i] = na*dotavx(a[1:p-i+1],a[i:p])
    end
    return b
end

function lpcacorr2v(rxx::AbstractVector{Float})
    k = 1/rxx[1]
    v = Vector{Float}(undef,length(rxx))
    for i ∈ eachindex(v)
        v[i] = rxx[i]*k
    end
    return v
end


function acorr_mat(rxx::AbstractVector{Float})
    p = length(rxx)
    V = Matrix{Float}(undef,p,p)
    for i ∈ 1:p
        for j ∈ 1:p
            V[j,i] = rxx[abs(i-j)+1]
        end
    end
    return V
end