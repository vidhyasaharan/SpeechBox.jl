#Running mean
function running_mean(x::AbstractMatrix{T}) where T<:AbstractFloat
    ndim = size(x,1)
    m = zeros(T,ndim)
    k = zero(T)
    for j ∈ axes(x,2)
        k = 1/j
        @turbo for i ∈ axes(x,1)
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
    @inbounds @views for i ∈ eachindex(mean)
        temp = (data[i] - mean[i])/n
        mean[i] += temp
    end
end