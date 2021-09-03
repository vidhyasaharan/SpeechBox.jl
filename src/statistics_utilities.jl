#Running mean
function running_mean(x::AbstractMatrix{T}) where T<:AbstractFloat
    ndim = size(x,1)
    m = zeros(T,ndim)
    k = zero(T)
    @inbounds for j ∈ axes(x,2)
        k = 1/j
        @inbounds @turbo for i ∈ axes(x,1)
            temp = (x[i,j] - m[i])
            m[i] += temp*k
        end
    end
    return m
end

#Running mean and variance
function running_meanvar(x::AbstractMatrix{T}) where T<:AbstractFloat
    ndim,npts = size(x)
    m = zeros(T,ndim)
    s = zeros(T,ndim)
    k = zero(T)
    for j ∈ axes(x,2)
        k = 1/j
        @turbo for i ∈ axes(x,1)
            temp = (x[i,j] - m[i])
            m[i] += temp*k
            s[i] += temp*(x[i,j] - m[i])
        end
    end
    N = 1/(npts-1)
    @turbo for i ∈ eachindex(s)
        s[i] *= N
    end
    return m, s
end


#Update running mean with a nth data point where n is the running index
function update_running_mean!(mean::AbstractVector{T}, data::AbstractVector{T}, n::Int) where T<:AbstractFloat
    k = 1/n
    @inbounds @turbo for i ∈ eachindex(mean)
        temp = (data[i] - mean[i])*k
        mean[i] += temp
    end
end


## Generic distance evaluation
#define Abstract type Distance
abstract type Distance end
struct SqL2 <: Distance end
struct L2 <: Distance end
dist(::SqL2, a::AbstractVector{T}, b::AbstractVector{T}) where {T} = sqL2avx(a,b)
dist(::L2, a::AbstractVector{T}, b::AbstractVector{T}) where {T} = L2avx(a,b)

#Pairwise distances
function pairwise!(dm::Distance, d::AbstractVector{<:AbstractFloat}, x::AbstractVector{T}, y::AbstractMatrix{T}) where {T}
    @inbounds @views for i ∈ axes(y,2)
        d[i] = dist(dm, x, y[:,i])
    end
end

function pairwise!(dm::Distance, d::AbstractMatrix{<:AbstractFloat}, x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T}
    @inbounds @views for i ∈ axes(x,2)
        for j ∈ axes(y,2)
            d[i,j] = dist(dm, x[:,i], y[:,j])
        end
    end
end

function pairwise(dm::Distance, x::AbstractVector{T}, y::AbstractMatrix{T}) where {T}
    d = Vector{Float}(undef,size(y,2))
    pairwise!(dm,d,x,y)
    return d
end

function pairwise(dm::Distance, x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T}
    d = Matrix{Float}(undef,size(x,2),size(y,2))
    pairwise!(dm,d,x,y)
    return d
end


#sample from a collection with a defined probability distribution

#get cumulative distribution (discrete) from scaled/unnormalised probability distribution (discrete)
function pdist2cdist(pdist::AbstractVector{T}) where {T}
    cdist = Vector{Float}(undef,length(pdist))
    cdist[1] = pdist[1]
    @views for i ∈ 2:length(pdist)
        cdist[i] = cdist[i-1] + pdist[i]
    end
    K = 1/cdist[end]
    @turbo for i ∈ eachindex(cdist)
        cdist[i] *= K
    end
    return cdist
end