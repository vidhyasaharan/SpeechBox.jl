##k-mean++

#Return distance to closest centre for each data point (helper function for k-mean++)
"""
    mindist2cntrs(centres, data)

Computes the square of the L2 distance to closest centre for each data point.
* `data` must be a matrix of `Float` with each column representing a point
* `centres` must be a matrix of `Float` with each column representing a centre

!!! note
    `mindist2cntrs(centres, data)` does not check the whether the dimensions of the input match.
    Users must ensure that `size(data,1) = size(centres,1)`.
"""
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
"""
    kmeanspp(data, number_of_centres)

Indentfies initial cluster centres for k-means clustering using the k-means++ alogrithm for given `data`.
* `data` must be a matrix of `Float` with each column vector representing a point
* Function will return a vector of indices corresponding to columns of `data` identified as cluster centres
"""
function kmeanspp(data::AbstractMatrix{T}, ncntrs::Int) where T<:AbstractFloat
    cin = [rand(1:size(data,2))]
    din = collect(1:size(data,2))
    deleteat!(din,cin)
    while(length(cin)<ncntrs)
        mindist = mindist2cntrs(data[:,cin], data[:,din])
        cdist = pdist2cdist(mindist)
        ncin = findclosest(rand(),cdist)
        push!(cin, din[ncin])
        deleteat!(din,ncin)
    end
    return cin
end


## K-means

#Identify closest cluster centre
"""
    closest_centre!(closest_centre, centres, data)

Identifies the centre closest to each data point based on L2 distance.
* `data` must be a matrix of `Float` with each column representing a point
* `centres` must be a matrix of `Float` with each column representing a centre
* `closest_centre` must be a vector of `Int` with each element holding the column index of `centres` that corresponds to the closest centre


!!! note
    `closest_centre!(closest_centre, centres, data)` does not check the whether the dimensions of the input match.
    Users must ensure that `length(closest_centre) = size(data,2)` and `size(data,1) = size(centres,1)`.
"""
function closest_centre!(ccntr::AbstractVector{Int}, cntrs::AbstractMatrix{T}, data::AbstractMatrix{T}) where T<:AbstractFloat
    temp = Vector{Float}(undef,size(cntrs,2))
    @inbounds @views for i ∈ eachindex(ccntr)
        pairwise!(SqL2(), temp, data[:,i], cntrs)
        ccntr[i] = argmin(temp)
    end
end

"""
    closest_centre(centres, data)

Identifies the centre closest to each data point based on L2 distance.
* `data` must be a matrix of `Float` with each column representing a point
* `centres` must be a matrix of `Float` with each column representing a centre
* Output is a vector of `Int` with each element holding the column index of `centres` that corresponds to the closest centre


!!! note
    `closest_centre(centres, data)` does not check the whether the dimensions of the input match.
    Users must ensure that `size(data,1) = size(centres,1)`.
"""
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

"""
    kmeans_init([initialise_method], data, number_of_centres)

Pick a set of points from `data` to serve as the initial cluster centres for k-means clustering. Method for initialisation can random or k-means++.
* `data` must be a matrix of `Float` with each column representing a point
* Use `kmpp()` as `initialise_method` to use k-mean++ as initialisation method [DEFAULT if `initialise_method` is not specified]
* Use `kmrand()` as `intialise_method` to randomly select points as initial centres. This can be much faster than k-means++ for large datasets
"""
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


"""
    kmeans!(centres, data, number_of_centres; [update_threshold = 0.0])

Updates (inplace) cluster centres `centres` using the k-means clustering algorithm based on `data`.
Cluster centres are iteratively updated until the largest shift in cluster centres,
measured in terms of square of L2 distance, falls below `update_threshold`.
* `centres` must be a matrix of `Float` with each column representing a cluster centre
* `data` must be a matrix of `Float` with each column representing a point
* Default value of optional named parameter is `update_threshold = 0.0` and corresponds to the setting that algorithm will converge when cluster centres do not shift
"""
function kmeans!(cntrs::AbstractMatrix{T}, data::AbstractMatrix{T}; update_threshold::Float = 0.0) where T<:AbstractFloat
    ndims,ncntrs = size(cntrs)
    old_cntrs = Matrix{T}(undef,ndims,ncntrs)
    ccntr = Vector{Int}(undef,size(data,2))
    max_shift::Float = Inf
    while(max_shift > update_threshold)
        copy!(old_cntrs,cntrs)
        kmeans_update!(ccntr, cntrs, data)
        max_shift = max_centre_shift(cntrs, old_cntrs)
    end
end


"""
    kmeans(data, number_of_centres; [update_threshold = 0.0])
    kmeans(data, number_of_centres, initialise_method; [update_threshold = 0.0])

Estimate cluster centres from `data` using the k-means clustering algorithm.
Cluster centres are iteratively updated until the largest shift in cluster centres,
measured in terms of square of L2 distance, falls below `update_threshold`.
* `data` must be a matrix of `Float` with each column representing a point
* Default value of optional named parameter is `update_threshold = 0.0` and corresponds to the setting that algorithm will converge when cluster centres do not shift
* Use `kmpp()` as `initialise_method` to use k-mean++ as initialisation method [DEFAULT if `initialise_method` is not specified]
* Use `kmrand()` as `intialise_method` to randomly select points as initial centres. This can be faster than k-means++ for large datasets
"""
function kmeans(data::AbstractMatrix{T}, ncntrs::Int, init_method::KMinit; update_threshold::Float = 0.0) where T<:AbstractFloat
    cntrs = kmeans_init(init_method, data, ncntrs)
    kmeans!(cntrs, data; update_threshold)
    return cntrs
end

kmeans(data::AbstractMatrix{T}, ncntrs::Int; update_threshold::Float = 0.0) where T<:AbstractFloat = kmeans(data, ncntrs, kmpp(); update_threshold)