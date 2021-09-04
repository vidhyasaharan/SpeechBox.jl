#Struct to hold parameters of a GMM (full covariance)
struct GMM
    w::Vector{Float}
    μ::Vector{Vector{Float}}
    Σ::Vector{Matrix{Float}}
end

## Constructors for GMM object

#Construct GMM when given a set of means, assuming equal weights and identity covariances
function GMM(means::AbstractMatrix{<:AbstractFloat})
    ndim,ncomp = size(means)
    w = convert(Vector{Float},(1/ncomp)*ones(ncomp))
    μ = [convert(Vector{Float},means[:,i]) for i ∈ axes(means,2)]
    Ī = convert(Matrix{Float},LinearAlgebra.I(ndim))
    Σ = [Ī for i ∈ 1:ncomp]
    return GMM(w,μ,Σ)
end


function SymmTri2Array(A::AbstractMatrix{T}) where {T}
    ndim = size(A,1)
    num_el::Int = ndim*(ndim+1)/2
    x = Vector{T}(undef,num_el)
    for i ∈ axes(A,1)
        x[i] = A[i,i]
    end
    n = ndim
    for i ∈ 1:ndim-1
        for j ∈ i+1:ndim
            n += 1
            x[n] = A[j,i] 
        end
    end
    return x
end
