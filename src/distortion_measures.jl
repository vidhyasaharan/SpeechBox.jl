
#Itakura Distortion
function distitak(x::AbstractVector{Float}, y::AbstractVector{Float}, p::Int)
    npts = 100
    Δθ = 2pi/(npts-1)
    θ = -pi:Δθ:pi
    pf1 = lpc_magz(x,2pi,p; frqs = θ)
    pf2 = lpc_magz(y,2pi,p; frqs = θ)
    d = zero(Float)
    for i ∈ eachindex(pf1)
        d += abs2(pf1[i]/pf2[i])
    end
    d *= Δθ/2pi
    return log(d)
end


function distitak2(x::AbstractVector{Float}, y::AbstractVector{Float}, p::Int)
    npts = 100
    Δθ = 2pi/(npts-1)
    θ = -pi:Δθ:pi
    rx = acorr(x,p+1)
    ry = acorr(y,p+1)
    αx,Ex,_ = levinson_durbin(rx)
    αy,Ey,_ = levinson_durbin(ry)
    ax = [1;-αx[end:-1:1]]
    ay = [1;-αy[end:-1:1]]

    fx = filter_coefs([1],ax)
    fy = filter_coefs([1],ay)

    hx = filter_magresp(fx, θ, 2pi)
    hy = filter_magresp(fy, θ, 2pi)

    pf1 = abs(Ex).*(hx.^2)
    pf2 = abs(Ey).*(hy.^2)

    λ = 0.01:0.01:2
    dis = Vector{Float}(undef,length(λ))
    for i ∈ eachindex(dis)
        dis[i] = distispf(pf1,λ[i].*pf2,Δθ)
    end
    return minimum(dis)
end


#Itakura Saito Distortion
function distis(x::AbstractVector{Float}, y::AbstractVector{Float}, p::Int)
    npts = 100
    Δθ = 2pi/(npts-1)
    θ = -pi:Δθ:pi
    pf1 = lpc_magz(x,2pi,p; frqs = θ)
    pf2 = lpc_magz(y,2pi,p; frqs = θ)
    d = distispf(pf1, pf2, Δθ)
    return d
end

function distispf(pf1::AbstractVector{Float}, pf2::AbstractVector{Float}, Δθ::Float)
    d = zero(Float)
    for i ∈ eachindex(pf1)
        v = pf1[i]/pf2[i]
        t = v - log(v) - 1
        d += t
    end
    d *= Δθ/2pi
    return d
end


# function distis(x::AbstractVector{Float}, y::AbstractVector{Float}, p::Int)
#     a = lpc(x, p)
#     b = lpcar2ra(a)
#     for i ∈ 2:length(b)
#         b[i] *= 2
#     end

#     Ryy = acorr(y, p+1)
#     â = lpc(Ryy)
#     r = lpcacorr2v(Ryy)

#     c = log(sum(abs2,a))
#     # d = c + log(dot(b,r)) - log(dot(â,r))
#     d = dot(b,r)/dot(â,r)
#     return d
# end


# function distis_mat(x::AbstractVector{Float}, y::AbstractVector{Float}, p::Int)
#     rxx = acorr(x, p+1)
#     V = acorr_mat(rxx)
#     a = lpc(x, p)
#     â = lpc(y, p)
#     N = dot(a, V*a)
#     D = dot(â, V*â)
#     d = N/D
#     return d
# end


#Convert inverse filter coefficients to autocorrelation coefficients
function lpcar2ra(a::AbstractVector{Float})
    na = 1/sum(abs2,a)
    p = length(a)
    b = Vector{Float}(undef,p)
    b[1] = 1
    for i ∈ 2:p
        b[i] = na*dot(a[1:p-i+1],a[i:p])
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