#Generate 2D data from 4 spherical (σ=0.25) gaussian clusters centered at [-1,1], [1,-1], [-1,1] and [1,1]
function generate_4_clusters(npts_per_cluster)
    c = convert(Matrix{Float},[1 1 -1 -1; 1 -1 1 -1])
    σ = 0.25
    data = Matrix{Float}(undef,2,4*npts_per_cluster)
    for i ∈ axes(c,2)
        rpts = randn(Float,(2,npts_per_cluster))
        sindx = (i-1)*npts_per_cluster+1
        eindx = i*npts_per_cluster
        data[:,sindx:eindx] = σ*rpts .+ c[:,i]
    end
    return data
end

##k-mean++

#Draw  Centre based on distances
# function kmpp_choose_cntr(min_dist::AbstractVector{T}) where {T}
#     r = rand()
# end

#Return distance to closest centre for each data point (helper function for k-mean++)
function mindist2cntrs(cntrs::AbstractMatrix{T}, data::AbstractMatrix{T}) where T<:AbstractFloat
    mdist = Vector{Float}(undef,size(data,2))
    temp = Vector{Float}(undef,size(cntrs,2))
    @inbounds @views for i ∈ eachindex(mdist)
        # temp = pairwise(SqL2(), data[:,i], cntrs)
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