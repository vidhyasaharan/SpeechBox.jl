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

# function update_running_mean(mean::AbstractVector{T}, data::AbstractVector{T}, data_index::Int) where T<:AbstractFloat
#     temp
# end