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
"""
    GMM(w::Vector{Float}, μ::Vector{Vector{Float}}, Σ::Vector{Matrix{Float}})

Create a `GMM` object given all its parameters.
* `w` must be a `Vector` of `Float` with the each element corresponding the weight of a Gaussian component
* `μ` must be a `Vector` of `Vector` of `Float`, each element of `μ` is a vector of `Float` representing the mean of a Gaussian component
* `Σ` must be a `Vector` of `Matrix` of `Float`, each element of `Σ` is a matrix of `Float` representing the covariance matrix of a Gaussian component

!!! note
    The elements of `w` must add up to 1.0 and the elements of `Σ` must be valid covariance matrices (symmetric and positive semi-definite) -  these
    constraints are not checked by the function. All elements of `μ` must be vectors of length `D` and all elements of `Σ` must be `D×D` matrices,
    where `D` is the dimensionality of the random variable being modelled.
    
"""
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
"""
    GMM(means)

Create a GMM object given the `means` of each component. The weights of all components are assumed to be equal and each covariance matrix is
assumed to be an identity matrix. The input `means` must be a matrix with each column representing the mean of a component.
    
"""
function GMM(means::AbstractMatrix{T}) where {T<:AbstractFloat}
    ndim,ncomp = size(means)
    w = convert(Vector{T},(1/ncomp)*ones(ncomp))
    μ = [means[:,i] for i ∈ axes(means,2)]
    Ī = convert(Matrix{T},LinearAlgebra.I(ndim))
    Σ = [Ī for i ∈ 1:ncomp]
    return GMM(w,μ,Σ)
end


#Initialise GMM parameters - equal weights, k-means for means and data covariance for each mixture covariance
abstract type GMinit end
struct init_kmeans <: GMinit end
struct init_rand <: GMinit end

GMMinit(nmix::Int, data::AbstractMatrix{T}) where {T<:AbstractFloat} = GMMinit(init_kmeans(), nmix, data)

function GMMinit(::init_kmeans, nmix::Int, data::AbstractMatrix{T}) where {T<:AbstractFloat}
    means = kmeans(data,nmix)
    m,C = running_meancov(data)
    w = ones(T,nmix)
    normaliseWeights!(w)
    μ = [means[:,i] for i ∈ axes(means,2)]
    Σ = [C for i ∈ 1:nmix]
    return GMM(w,μ,Σ)
end

function GMMinit(::init_rand, nmix::Int, data::AbstractMatrix{T}) where {T<:AbstractFloat}
    ndims,npts = size(data)
    mindx = rand(1:npts,nmix)
    means = data[:,mindx]
    II = convert(Matrix{Float},collect(I(ndims)))
    w = ones(T,nmix)
    normaliseWeights!(w)
    μ = [means[:,i] for i ∈ axes(means,2)]
    Σ = [II for i ∈ 1:nmix]
    return GMM(w,μ,Σ)
end

#Compute Log Probability (Log-likelihood)
logprob(μ::AbstractVector{T}, P::AbstractMatrix{T}, Z::T, x::AbstractVector{T}) where {T<:AbstractFloat} = Z-(mahalavx(x,μ,P)/2)

function logmixprob!(lprobs::AbstractVector{Float}, G::GMM, x::AbstractVector{Float})
    for i ∈ eachindex(lprobs, G.μ, G.P, G.Z, G.w)
        lprobs[i] = logprob(G.μ[i], G.P[i], G.Z[i], x) + log(G.w[i])
    end
end

"""
    logprob(G::GMM, x)

Compute the log probability of data `x` given a Gaussian mixture model `G` (provided as a `GMM` object). Equivalently,
this is the log-likelyhood of model `G` given the data `x`.
* when `x` is a vector representing one data point, the log probability of the data point given the model is returned
* when `x` is a `D×N` matrix representing a set of `N` points of dimensionality `D`, the sum of the log probabilities
  of each point is returned (i.e., total probability assuming points are independent)
"""
function logprob(G::GMM, x::AbstractVector{Float})
    nmix = length(G.w)
    lprobs = Vector{Float}(undef,nmix)
    logmixprob!(lprobs,G,x)
    return logsumexp(lprobs)
end

function logprob(G::GMM, x::AbstractMatrix{Float})
    LL = zero(Float)
    @views for i ∈ axes(x,2)
        LL += logprob(G,x[:,i])
    end
    return LL
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
"""
    normaliseWeights!(weights)

Normalises the elements of the vector `weights` in place such that they sum to 1.0
"""
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

#Sample points from a GMM
"""
    sample(G::GMM, N::Integer)

Draw `N` sample from a Gaussian mixture model `G` provided as a `GMM` object
"""
function sample(G::GMM, npts::Int)
    data = Matrix{Float}(undef,length(G.μ[1]),npts)
    gcomps = Vector{Gaussian}(undef,length(G.w))
    wdist = Categorical(G.w)
    for i ∈ eachindex(G.μ,G.Σ)
        gcomps[i] = Gaussian(G.μ[i], G.Σ[i])
    end
    for i ∈ axes(data,2)
        cin = sample(wdist)
        data[:,i] = sample(gcomps[cin])
    end
    return data
end


#Struct to hold categorical distribution
struct Categorical
    pdist::Vector{Float}
    cdist::Vector{Float}
end

"""
    Categorical(probability_vector)

Create a `Categorical` object to hold probability and cumulative distributions of a finite set of discrete outcomes
given `probability_vector` a `Vector` of probabilities of each outcome.

!!! note
    The elements of `probability_vector` should sum to 1.0. This is NOT checked.
"""
function Categorical(pdist::AbstractVector{Float})
    cdist = pdist2cdist(pdist)
    return Categorical(pdist,cdist)
end

"""
    sample(C::Categorical)

Draw a sample from a categorical distribution `C` provided as a `Categorical` object
"""
sample(C::Categorical) = findclosest(rand(), C.cdist)


#Struct to hold a multivariate Gaussian
struct Gaussian
    μ::Vector{Float}
    Σ::Matrix{Float}
    A::Cholesky{Float}
end

"""
    Gaussian(μ, Σ)

Create a `Gaussian` object to hold the parameters of a Gaussian distribution when its mean `μ` is given as a vector (of Floats) and its
covariance matrix `Σ` is given as a matrix (of Floats)

!!! note
    `Σ` must be a symmetric, positive semi-definite matrix. This is NOT checked.
"""
Gaussian(μ::AbstractVector{T}, Σ::AbstractMatrix{T}) where {T<:AbstractFloat} = Gaussian(μ,Σ,cholesky(Σ))

"""
    sample(g::Gaussian)
    sample(g::Gaussian, N::Integer)

Draw one or `N` samples from a Gaussian distributions `g` provided as a `Gaussian` object
"""
function sample(g::Gaussian)
    z = randn(Float,length(g.μ))
    return g.μ + g.A.L*z
end

function sample(g::Gaussian, npts::Int)
    data = Matrix{Float}(undef,length(g.μ),npts)
    for i ∈ axes(data,2)
        data[:,i] = sample(g)
    end
    return data
end

function generate_4mix_GMM()
    c = 4*convert(Matrix{Float},[1 1 -1 -1; 1 -1 1 -1])
    μ = [c[:,i] for i ∈ axes(c,2)]
    Σ₁ = [1 0; 0 1]
    Σ₂ = [2 0; 0 1]
    Σ₃ = [1 .75; .75 1]
    Σ₄ = [1 -.75; -.75 1]
    Σ = convert(Vector{Matrix{Float}},[Σ₁, Σ₂, Σ₃, Σ₄])
    w = [0.1, 0.2, 0.3, 0.4]
    return GMM(w,μ,Σ)
end

## EM Algorithm (ML Estimate)

#E-step
function mixture_posterior!(γ::AbstractVector{Float}, G::GMM, x::AbstractVector{Float})
    logmixprob!(γ,G,x)
    lp = logsumexp(γ)
    @turbo for i ∈ eachindex(γ)
        γ[i] = exp(γ[i] - lp)
    end
end

function mixture_posterior(G::GMM, x::AbstractVector{Float})
    γ = Vector{Float}(undef,length(G.w))
    mixture_posterior!(γ,G,x)
    return γ
end

function mixture_posterior!(γ::AbstractMatrix{Float}, G::GMM, x::AbstractMatrix{Float})
    @views for i ∈ axes(x,2)
        mixture_posterior!(γ[:,i],G,x[:,i])
    end
end

function mixture_posterior(G::GMM, x::AbstractMatrix{Float})
    γ = Matrix{Float}(undef,length(G.w), size(x,2))
    mixture_posterior!(γ,G,x)
    return γ
end

#M-step
function update_ML(G::GMM, γ::AbstractMatrix{Float}, x::AbstractMatrix{Float})
    N = size(x,2)
    ndim = length(G.μ[1])
    nmix = length(G.w)
    w = Vector{Float}(undef,nmix)
    μ = [zeros(Float,ndim) for i ∈ 1:nmix]
    Σ = [zeros(Float,ndim,ndim) for i ∈ 1:nmix]
    k = Vector{Float}(undef,nmix)
    temp = Matrix{Float}(undef,ndim,nmix)
    Nm = zeros(Float,nmix)
    for j ∈ axes(x,2)
        for i ∈ 1:nmix
            Nm[i] += γ[i,j]
        end
    end
    for i ∈ eachindex(w,Nm)
        w[i] = Nm[i]/N
        k[i] = 1/Nm[i]
    end
    for m ∈ 1:nmix
        for n ∈ axes(x,2)
            kn = γ[m,n]*k[m]
            @turbo for i ∈ 1:ndim
                μ[m][i] += kn*x[i,n]
                temp[i,m] = x[i,n] - G.μ[m][i]
            end
            @turbo for i ∈ 1:ndim, j ∈ 1:ndim
                Σ[m][i,j] += kn*temp[i,m]*temp[j,m]
            end
        end
        Σ[m] += Σ[m]'
        Σ[m] /= 2
    end
    return GMM(w,μ,Σ)
end

"""
    trainML(G::GMM, data, num_iterations)
    trainML([initalise_method], data, num_iterations)

Train a Maximum Likelihood (ML) Gaussian mixture model fit for `data` using the Expectation-Maximisation (EM) algorithm.
`data` must be a matrix of `Float` with each column representing a point. Initial values for the Gaussian mixture model
parameters may be obtained from a GMM object `G` or estimated from `data` using a suitable initialisation method:
* Use `init_kmeans()` as `initialise_method` to use k-means (with k-mean++ initialisation) on `data` to initialise the means, the
  weights are set as equal weights, covariance matrices of each component are set to covariance of `data`
* Use `init_rand()` as `initialise_method` to select poitsn from `data` at random as means, the weights are set as equal weights, and
  the covariance matrices of each component are set as identity matrices
"""
function trainML(G::GMM, x::AbstractMatrix{Float}, niter::Int)
    γ = Matrix{Float}(undef,length(G.w),size(x,2))
    LL = Vector{Float}(undef,0)
    for i ∈ 1:niter
        mixture_posterior!(γ,G,x)
        G = update_ML(G, γ, x)
        push!(LL,logprob(G,x))
    end
    return G, LL
end

function trainML(itype::GMinit, x::AbstractMatrix{Float}, nmix::Int, niter::Int)
    G = GMMinit(itype, nmix, x)
    return trainML(G, x, niter)
end

trainML(x::AbstractMatrix{Float}, nmix::Int, niter::Int) = trainML(init_kmeans(), x, nmix, niter)