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
    fill!(m, zero(T))
    @inbounds for j ∈ axes(x,2)
        k = 1/j
        @simd for i ∈ axes(x,1)
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
    fill!(m, zero(T))
    fill!(v, zero(T))
    @inbounds for j ∈ axes(x,2)
        k = 1/j
        for i ∈ axes(x,1)
            temp = (x[i,j] - m[i])
            m[i] += temp*k
            v[i] += temp*(x[i,j] - m[i])
        end
    end
    N = 1/(size(x,2)-1)
    v .*= N
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
        @simd for i ∈ axes(x,1)
            temp[i] = (x[i,n] - m[i])
            m[i] += temp[i]*k
        end
        for i ∈ axes(x,1)
            @simd for j ∈ axes(x,1)
                C[j,i] += (x[j,n] - m[j])*temp[i]
            end
        end
    end
    N = 1/(size(x,2)-1)
    C .*= N
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
    @inbounds @simd for i ∈ eachindex(mean)
        temp = (data[i] - mean[i])*k
        mean[i] += temp
    end
end


## Distance computations

#Squared L2 (Euclidean) distance between two vectors
function sqL2dist(a::AbstractVector{T}, b::AbstractVector{T}) where {T}
    s = zero(T)
    @inbounds @simd for i ∈ eachindex(a,b)
        t = a[i] - b[i]
        s += t * t
    end
    return s
end

#L2 (Euclidean) distance between two vectors
L2dist(a::AbstractVector{T}, b::AbstractVector{T}) where {T} = sqrt(sqL2dist(a,b))

#Squared Mahalanobis-type distance (x-y)'A(x-y), where A is typically a precision matrix
function sqmahal(x::AbstractVector{T}, y::AbstractVector{T}, A::AbstractMatrix{T}) where {T}
    s = zero(T)
    @inbounds for j ∈ eachindex(x,y)
        tj = x[j] - y[j]
        r = zero(T)
        @simd for i ∈ eachindex(x,y)
            r += (x[i] - y[i]) * A[i,j]
        end
        s += tj * r
    end
    return s
end


## Generic distance evaluation
#define Abstract type Distance
abstract type Distance end
struct SqL2 <: Distance end
struct L2 <: Distance end
dist(::SqL2, a::AbstractVector{T}, b::AbstractVector{T}) where {T} = sqL2dist(a,b)
dist(::L2, a::AbstractVector{T}, b::AbstractVector{T}) where {T} = L2dist(a,b)

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
    cdist .*= K
    return cdist
end