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

##k-mean++

#Return distance to closest centre for each data point (helper function for k-mean++)
function mindist2cntrs(cntrs::AbstractMatrix{T}, data::AbstractMatrix{T}) where T<:AbstractFloat
    mdist = Vector{Float}(undef,size(data,2))
    temp = Vector{Float}(undef,size(cntrs,2))
    @inbounds @views for i ∈ eachindex(mdist)
        pairwise!(SqL2(), temp, data[:,i], cntrs)
        mdist[i] = minimum(temp)
    end
    return mdist
end

#k-means++ algorithm for initialising centres
function kmeanspp(data::AbstractMatrix{T}, ncntrs::Int) where T<:AbstractFloat
    cin = [rand(1:size(data,2))]
    din = collect(1:size(data,2))
    deleteat!(din,cin)
    while(length(cin)<ncntrs)
        mindist = SpeechBox.mindist2cntrs(data[:,cin], data[:,din])
        cdist = SpeechBox.pdist2cdist(mindist)
        ncin = SpeechBox.findclosest(rand(),cdist)
        push!(cin, din[ncin])
        deleteat!(din,ncin)
    end
    return cin
end


## K-means

#Identify closest cluster centre
function closest_centre!(ccntr::AbstractVector{Int}, cntrs::AbstractMatrix{T}, data::AbstractMatrix{T}) where T<:AbstractFloat
    temp = Vector{Float}(undef,size(cntrs,2))
    @inbounds @views for i ∈ eachindex(ccntr)
        pairwise!(SqL2(), temp, data[:,i], cntrs)
        ccntr[i] = argmin(temp)
    end
end

function closest_centre(cntrs::AbstractMatrix{T}, data::AbstractMatrix{T}) where T<:AbstractFloat
    ccntr = Vector{Int}(undef,size(data,2))
    closest_centre!(ccntr, cntrs, data)
    return ccntr
end

#Maximum L2 distance by which a centre has shifted
function max_centre_shift(new_cntrs::AbstractMatrix{T}, old_cntrs::AbstractMatrix{T}) where T<:AbstractFloat
    s = zero(T)
    @views for i ∈ axes(old_cntrs,2)
        temp = dist(SqL2(),new_cntrs[:,i],old_cntrs[:,i])
        if(temp>s)
            s = temp
        end
    end
    return sqrt(s)
end

#k-means initialisation
abstract type KMinit end
struct kmpp <: KMinit end
struct kmrand <: KMinit end

kmeans_init(::kmpp, data::AbstractMatrix{<:AbstractFloat}, ncntrs::Int) = data[:,kmeanspp(data,ncntrs)]

kmeans_init(::kmrand, data::AbstractMatrix{<:AbstractFloat}, ncntrs::Int) = data[:,rand(1:size(data,2),ncntrs)]

kmeans_init(data::AbstractMatrix{<:AbstractFloat}, ncntrs::Int) = kmeans_init(kmpp(), data, ncntrs)


#One iteration of updating k-means cluster centres
function kmeans_update!(ccntr::AbstractVector{Int}, cntrs::AbstractMatrix{T}, data::AbstractMatrix{T}) where T<:AbstractFloat
    closest_centre!(ccntr,cntrs,data)
    for i ∈ axes(cntrs,2)
        if(sum(ccntr.==i)>0)
            cntrs[:,i] = running_mean(data[:,ccntr.==i])
        end
    end
end



function kmeans!(cntrs::AbstractMatrix{T}, data::AbstractMatrix{T}; threshold::Float = 0.0) where T<:AbstractFloat
    ndims,ncntrs = size(cntrs)
    old_cntrs = Matrix{T}(undef,ndims,ncntrs)
    ccntr = Vector{Int}(undef,size(data,2))
    max_shift::Float = Inf
    while(max_shift > threshold)
        copy!(old_cntrs,cntrs)
        kmeans_update!(ccntr, cntrs, data)
        max_shift = max_centre_shift(cntrs, old_cntrs)
    end
end

function kmeans(data::AbstractMatrix{T}, ncntrs::Int; threshold::Float = 0.0) where T<:AbstractFloat
    cntrs = kmeans_init(kmpp(), data, ncntrs)
    kmeans!(cntrs, data; threshold)
    return cntrs
end