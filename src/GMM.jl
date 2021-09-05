#Struct to hold parameters of a GMM (full covariance)
struct GMM
    w::Vector{Float}
    μ::Vector{Vector{Float}}
    Σ::Vector{Matrix{Float}}
    P::Vector{Matrix{Float}} #Precision matrix P = inv(Σ)
    Z::Vector{Float} #Normalising constant for log probability of Gaussian: Z = -(D/2)log(2π)-(1/2)log|Σ|
end

## Constructors for GMM object

#Construct GMM given a set of weights, means and covariance matrices
function GMM(w::AbstractVector{T}, μ::Vector{Vector{T}}, Σ::Vector{Matrix{T}}) where {T<:AbstractFloat}
    if(isconsistentParams(w,μ,Σ))
        D::Int = length(μ[1])
        Z₀::Float = -(D/2)*log(2π)
        P = [inv(C) for C ∈ Σ]
        Z = Vector{Float}(undef,length(w))
        for i ∈ eachindex(Z,Σ)
            Z[i] = Z₀ - (logdet(Σ[i])/2)
        end
        return GMM(w,μ,Σ,P,Z)
    else
        return nothing
    end
end

#Construct GMM when given a set of means, assuming equal weights and identity covariances
function GMM(means::AbstractMatrix{T}) where {T<:AbstractFloat}
    ndim,ncomp = size(means)
    w = convert(Vector{T},(1/ncomp)*ones(ncomp))
    μ = [means[:,i] for i ∈ axes(means,2)]
    Ī = convert(Matrix{T},LinearAlgebra.I(ndim))
    Σ = [Ī for i ∈ 1:ncomp]
    return GMM(w,μ,Σ)
end


#Initialise GMM parameters - equal weights, k-means for means and data covariance for each mixture covariance
function GMMinit(nmix::Int, data::AbstractMatrix{T}) where {T<:AbstractFloat}
    means = kmeans(data,nmix)
    m,C = running_meancov(data)
    w = ones(T,nmix)
    normaliseWeights!(w)
    μ = [means[:,i] for i ∈ axes(means,2)]
    Σ = [C for i ∈ 1:nmix]
    return GMM(w,μ,Σ)
end

#Compute Log Probability (Log-likelihood)
logprob(μ::AbstractVector{T}, P::AbstractMatrix{T}, Z::T, x::AbstractVector{T}) where {T<:AbstractFloat} = Z-(mahalavx(x,μ,P)/2)

function logprob(G::GMM, x::AbstractVector{Float})
    nmix = length(G.w)
    lprobs = Vector{Float}(undef,nmix)
    for i ∈ eachindex(lprobs, G.μ, G.P, G.Z, G.w)
        lprobs[i] = logprob(G.μ[i], G.P[i], G.Z[i], x) + log(G.w[i])
    end
    return logsumexp(lprobs)
end

function logsumexp(lp::AbstractVector{T}) where T<:AbstractFloat
    lmax = maximum(lp)
    s = zero(T)
    @turbo for i ∈ eachindex(lp)
        s += exp(lp[i] - lmax)
    end
    return log(s) + lmax
end


## Utility functions to check consistency among GMM parameters

#Normalise component weights so that they sum to one
function normaliseWeights!(w::AbstractVector{<:AbstractFloat})
    w₀ = 1/sum(w)
    for i ∈ eachindex(w)
        w[i] *= w₀
    end
end

#Combine all consistency checks in one function
function isconsistentParams(w::AbstractVector{T}, μ::Vector{Vector{T}}, Σ::Vector{Matrix{T}}) where {T<:AbstractFloat}
    if(!isconsistentComponents(w,μ,Σ))
        return false
    elseif(!isconsistentWeights(w))
        return false
    elseif(!isconsistentDimensions(μ,Σ))
        return false
    end
    return true
end

#Check if Number of components are consistent
function isconsistentComponents(w::AbstractVector{T}, μ::Vector{Vector{T}}, Σ::Vector{Matrix{T}}) where {T<:AbstractFloat}
    nc_w = length(w)
    nc_μ = length(μ)
    nc_Σ = length(Σ)
    wμ = (nc_w == nc_μ)
    wΣ = (nc_w == nc_Σ)
    μΣ = (nc_μ == nc_Σ)
    if(wμ & wΣ & μΣ)
        return true
    else
        return false
    end
end

#Check if weights sum to one
function isconsistentWeights(w::AbstractVector{T}) where {T<:AbstractFloat}
    if(sum(w) ≈ one(T))
        return true
    else
        return false
    end
end

#Check if dimensions of all means and covariances are consistent
function isconsistentDimensions(μ::Vector{Vector{T}}, Σ::Vector{Matrix{T}}) where {T<:AbstractFloat}
    D = length(μ[1])
    flag::Bool = true
    for i ∈ eachindex(μ,Σ)
        if(isempty(μ[i])||isempty(Σ[i]))
            return false
        end
        length(μ[i])==D ? flag &= true : flag &= false
        size(Σ[i],1)==D ? flag &= true : flag &= false
        size(Σ[i],2)==D ? flag &= true : flag &= false
    end
    return flag
end


