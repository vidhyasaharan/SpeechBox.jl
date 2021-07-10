
function acorr!(rxx::AbstractVector{Float}, x::AbstractVector{Float})
    len = length(x)
    @inbounds @views for i ∈ eachindex(rxx)
        δ = i-1
        N = len - δ
        rxx[i] = (1/N)*dotavx(x[1:end-δ],x[1+δ:end])
    end
end

function acorr(x::AbstractVector{Float}, p::Int)
    rxx = Vector{Float}(undef,p)
    acorr!(rxx,x)
    return rxx
end


function lpc_LD(x::AbstractVector{Float}, p::Int)
    α = Vector{Float}(undef,p)
    lpc_LD!(α,x)
    return [1;-α[end:-1:1]]
end


function lpc_LD!(α::AbstractVector{Float}, x::AbstractVector{Float})
    p = length(α)
    rxx = acorr(x, p+1)
    E = Vector{Float}(undef,p+1)
    k = Vector{Float}(undef,p)
    @views E[1] = rxx[1]
    @views k[1] = rxx[2]/rxx[1]
    @views α[p] = k[1]
    update_E!(E,k,1)
    for i ∈ 2:p
        update_k!(k, α, rxx, E, i)
        update_α!(α, k, i)
        update_E!(E, k, i)
    end
end

function update_E!(E::AbstractVector{Float}, k::AbstractVector{Float}, iter::Int)
    @views E[iter+1] = (1-(k[iter]^2))*E[iter]
end

function update_k!(k::AbstractVector{Float}, α::AbstractVector{Float}, rxx::AbstractVector{Float}, E::AbstractVector{Float}, iter::Int)
    @views a = α[end-iter+2:end]
    @views k[iter] = (rxx[iter+1] - dotavx(a,rxx[2:iter]))/E[iter]
end

function update_α!(α::AbstractVector{Float}, k::AbstractVector{Float}, iter::Int)
    @views α[end-iter+1] = k[iter]
    a = α[end-iter+2:end]
    @views for i ∈ 1:iter-1
        α[end-i+1] -= k[iter]*a[i]
    end
end
