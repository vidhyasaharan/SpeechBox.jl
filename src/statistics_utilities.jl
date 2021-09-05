#Generate 2D data from 4 spherical (σ=0.25) gaussian clusters centered at [-1,1], [1,-1], [-1,1] and [1,1]
function generate_4_clusters(npts_per_cluster; σ = 0.25)
    c = convert(Matrix{Float},[1 1 -1 -1; 1 -1 1 -1])
    # σ = 0.25
    data = Matrix{Float}(undef,2,4*npts_per_cluster)
    for i ∈ axes(c,2)
        rpts = randn(Float,(2,npts_per_cluster))
        sindx = (i-1)*npts_per_cluster+1
        eindx = i*npts_per_cluster
        data[:,sindx:eindx] = σ*rpts .+ c[:,i]
    end
    return data
end

#Generate 2D data arranged in a circle (default radius, r =0.5) around centres at [-1,1], [1,-1], [-1,1] and [1,1]
function generate_4_circle_clusters(npts_per_cluster::Int = 8; r::Float = 0.5)
    c = convert(Matrix{Float},[1 1 -1 -1; 1 -1 1 -1])
    data = Matrix{Float}(undef,2,4*npts_per_cluster)
    for i ∈ axes(c,2)
        sindx = (i-1)*npts_per_cluster+1
        eindx = i*npts_per_cluster
        pts = generate_circle_cluster(npts_per_cluster)
        data[:,sindx:eindx] = r*pts .+ c[:,i]
    end
    return data
end

#Generate 2D data arranged in a circle of radius one around [0,0]
function generate_circle_cluster(npts::Int = 8)
    Δθ = 2π/npts
    θ = 0:Δθ:2π-Δθ
    x = Matrix{Float}(undef,2,npts)
    for i ∈ eachindex(θ)
        x[1,i] = cos(θ[i])
        x[2,i] = sin(θ[i])
    end
    return x
end

#Running mean
function running_mean!(m::AbstractVector{T}, x::AbstractMatrix{T}) where T<:AbstractFloat
    k = zero(T)
    @turbo for i ∈ eachindex(m)
        m[i] = zero(T)
    end
    @inbounds for j ∈ axes(x,2)
        k = 1/j
        @turbo for i ∈ axes(x,1)
            temp = (x[i,j] - m[i])
            m[i] += temp*k
        end
    end
end

function running_mean(x::AbstractMatrix{T}) where T<:AbstractFloat
    m = zeros(T,size(x,1))
    running_mean!(m, x)
    return m
end

#Running mean and variance
function running_meanvar!(m::AbstractVector{T}, v::AbstractVector{T}, x::AbstractMatrix{T}) where T<:AbstractFloat
    k = zero(T)
    @turbo for i ∈ eachindex(m)
        m[i] = zero(T)
        v[i] = zero(T)
    end
    @inbounds for j ∈ axes(x,2)
        k = 1/j
        @turbo for i ∈ axes(x,1)
            temp = (x[i,j] - m[i])
            m[i] += temp*k
            v[i] += temp*(x[i,j] - m[i])
        end
    end
    N = 1/(size(x,2)-1)
    @turbo for i ∈ eachindex(v)
        v[i] *= N
    end
end

function running_meanvar(x::AbstractMatrix{T}) where T<:AbstractFloat
    m = zeros(T,size(x,1))
    v = zeros(T,size(x,1))
    running_meanvar!(m,v,x)
    return m, v
end


#running mean and covariance
function running_meancov!(m::AbstractVector{T}, C::AbstractMatrix{T}, x::AbstractMatrix{T}) where T<:AbstractFloat
    k = zero(T)
    temp = Vector{T}(undef,length(m))
    fill!(m,zero(T))
    fill!(C,zero(T))
    @inbounds for n ∈ axes(x,2)
        k = 1/n
        @turbo for i ∈ axes(x,1)
            temp[i] = (x[i,n] - m[i])
            m[i] += temp[i]*k
        end
        @turbo for i ∈ axes(x,1), j ∈ axes(x,1)
            C[j,i] += (x[j,n] - m[j])*temp[i]
        end
    end
    N = 1/(size(x,2)-1)
    @turbo for i ∈ eachindex(C)
        C[i] *= N
    end
    C += C'
    C /= 2
end

function running_meancov(x::AbstractMatrix{T}) where T<:AbstractFloat
    m = Vector{Float}(undef,size(x,1))
    C = Matrix{Float}(undef,size(x,1),size(x,1))
    running_meancov!(m,C,x)
    return m, C
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