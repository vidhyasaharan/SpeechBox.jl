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

function logmixprob!(lprobs::AbstractVector{Float}, G::GMM, x::AbstractVector{Float})
    for i ∈ eachindex(lprobs, G.μ, G.P, G.Z, G.w)
        lprobs[i] = logprob(G.μ[i], G.P[i], G.Z[i], x) + log(G.w[i])
    end
end

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

function Categorical(pdist::AbstractVector{Float})
    cdist = pdist2cdist(pdist)
    return Categorical(pdist,cdist)
end

sample(C::Categorical) = findclosest(rand(), C.cdist)


#Struct to hold a multivariate Gaussian
struct Gaussian
    μ::Vector{Float}
    Σ::Matrix{Float}
    A::Cholesky{Float}
end

Gaussian(μ::AbstractVector{T}, Σ::AbstractMatrix{T}) where {T<:AbstractFloat} = Gaussian(μ,Σ,cholesky(Σ))

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