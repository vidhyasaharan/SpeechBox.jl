function posdefmatrix(ndim::Int)
    a = randn(Float64,ndim,ndim)
    A = a*a' + LinearAlgebra.I(ndim)
    return A
end

function normpdf(x::AbstractVector{T}, μ::AbstractVector{T}, Σ::AbstractMatrix{T}) where {T<:AbstractFloat}
    x̄ = x - μ
    D = length(x)
    z = sqrt(det(Σ)*((2π)^D))
    return (exp(-(x̄'*inv(Σ)*x̄)/2))/z
end

@testset "normaliseWts" begin
    nmix = 10
    w = rand(Float64,nmix)
    ŵ = copy(w)
    SpeechBox.normaliseWeights!(w)
    @test sum(w) ≈ one(Float64)
    for i ∈ eachindex(w,ŵ)
        @test w[i]/sum(w) ≈ ŵ[i]/sum(ŵ)
    end
end

@testset "isconsistentX" begin
    nmix = 5
    ndim = 3
    w = rand(nmix)
    SpeechBox.normaliseWeights!(w)
    μ = [randn(Float64,ndim) for i ∈ 1:nmix]
    Σ = [posdefmatrix(ndim) for i ∈ 1:nmix]

    @test SpeechBox.isconsistentComponents(w,μ,Σ) == true
    @test SpeechBox.isconsistentComponents(w[1:nmix-1],μ,Σ) == false
    @test SpeechBox.isconsistentComponents(w,μ[1:nmix-1],Σ) == false
    @test SpeechBox.isconsistentComponents(w,μ,Σ[1:nmix-1]) == false

    @test SpeechBox.isconsistentWeights(w) == true
    @test SpeechBox.isconsistentWeights(ones(nmix)) == false

    @test SpeechBox.isconsistentDimensions(μ,Σ) == true
    μ̄ = copy(μ)
    μ̄[nmix] = randn(Float64,ndim-1)
    Σ̄ = copy(Σ)
    Σ̄[1] = posdefmatrix(ndim+1)
    @test SpeechBox.isconsistentDimensions(μ̄,Σ) == false
    @test SpeechBox.isconsistentDimensions(μ,Σ̄) == false

    @test SpeechBox.isconsistentParams(w,μ,Σ) == true
    @test SpeechBox.isconsistentParams(w[1:nmix-1],μ,Σ) == false
    @test SpeechBox.isconsistentParams(w,μ[1:nmix-1],Σ) == false
    @test SpeechBox.isconsistentParams(w,μ,Σ[1:nmix-1]) == false
    @test SpeechBox.isconsistentParams(ones(nmix),μ,Σ) == false
    @test SpeechBox.isconsistentParams(w,μ̄,Σ) == false
    @test SpeechBox.isconsistentParams(w,μ,Σ̄) == false
end


@testset "GMM Struct" begin
    ndim = 10
    nmix = 5
    w = rand(nmix)
    SpeechBox.normaliseWeights!(w)
    μ = [randn(Float64,ndim) for i ∈ 1:nmix]
    Σ = [posdefmatrix(ndim) for i ∈ 1:nmix]

    G = SpeechBox.GMM(w,μ,Σ)

    @test typeof(G.w) == Vector{Float64}
    @test typeof(G.μ) == Vector{Vector{Float64}}
    @test typeof(G.Σ) == Vector{Matrix{Float64}}
    @test typeof(G.P) == Vector{Matrix{Float64}}
    @test typeof(G.Z) == Vector{Float64}

    @test length(G.w) == nmix
    @test G.w == w

    
    for i = 1:nmix
        @test G.μ[i] == μ[i]
        @test G.Σ[i] == Σ[i]
        @test G.P[i] == inv(Σ[i])
        @test G.Z[i] + (ndim/2)*log(2π) ≈ -logdet(Σ[i])/2
    end


    means = randn(ndim,nmix)
    G = SpeechBox.GMM(means)

    @test typeof(G.w) == Vector{Float64}
    @test typeof(G.μ) == Vector{Vector{Float64}}
    @test typeof(G.Σ) == Vector{Matrix{Float64}}
    @test typeof(G.P) == Vector{Matrix{Float64}}
    @test typeof(G.Z) == Vector{Float64}

    @test G.w == (1/nmix)*ones(nmix)
    for i = 1:nmix
        @test G.μ[i] == means[:,i]
        @test G.Σ[i] == I(ndim)
        @test G.P[i] == I(ndim)
    end
end

@testset "GMMinit" begin
    x = SpeechBox.generate_4_circle_clusters(r=0.1)
    kmeans = SpeechBox.kmeans(x,4)
    G = SpeechBox.GMMinit(SpeechBox.init_kmeans(),4,x)
    @test G.w == (1/4)*ones(4)
    @test SpeechBox.mindist2cntrs(kmeans,hcat(G.μ...)) == zeros(4)

    G = SpeechBox.GMMinit(SpeechBox.init_rand(),4,x)
    @test SpeechBox.mindist2cntrs(x, hcat(G.μ...)) == zeros(4)
end

@testset "logsumexp" begin
    nprbs = 10
    probs = rand(Float64,nprbs)
    lprobs = log.(probs)
    @test SpeechBox.logsumexp(lprobs) ≈ log(sum(probs))
end

@testset "Categorical" begin
    ncat = 10
    p = rand(Float64,ncat)
    SpeechBox.normaliseWeights!(p)
    # c = SpeechBox.pdist2cdist(p)
    c = SpeechBox.Categorical(p)
    @test c isa SpeechBox.Categorical{Float64}
    for i ∈ eachindex(p,c.pdist,c.cdist)
        @test c.pdist[i] ≈ p[i]
        @test c.cdist[i] ≈ sum(p[1:i])
    end
end

@testset "Gaussian" begin
    ndim = 5
    μ = randn(Float64,ndim)
    Σ = posdefmatrix(ndim)
    g = SpeechBox.Gaussian(μ,Σ)
    @test g isa SpeechBox.Gaussian{Float64}
    @test μ == g.μ
    @test Σ == g.Σ
    @test LinearAlgebra.cholesky(Σ) == g.A
    @test g.A.L*g.A.U ≈ g.Σ
end


@testset "logprob" begin
    D = 5
    μ = randn(D)
    Σ = posdefmatrix(D)
    P = inv(Σ)
    Z = -(D/2)*log(2π) - (logdet(Σ)/2)
    for i ∈ 1:10
        x = randn(D)
        @test SpeechBox.logprob(μ,P,Z,x) ≈ log(normpdf(x,μ,Σ))
    end

    G = SpeechBox.generate_4mix_GMM()
    for i ∈ 1:10
        x = randn(2)
        p = 0
        for j ∈ 1:4
            p += G.w[j]*normpdf(x,G.μ[j],G.Σ[j])
        end
        @test SpeechBox.logprob(G,x) ≈ log(p)
    end

    for n ∈ 1:10
        x = randn(2,10)
        ll = 0
        for i ∈ axes(x,2)
            ll += SpeechBox.logprob(G,x[:,i])
        end
        @test SpeechBox.logprob(G,x) ≈ ll
    end
end

@testset "mix_posterior" begin
    G = SpeechBox.generate_4mix_GMM()
    for m ∈ 1:4
        x = SpeechBox.sample(SpeechBox.Gaussian(G.μ[m],G.Σ[m]),10)
        γ = SpeechBox.mixture_posterior(G,x)
        for i ∈ axes(γ,2)
            @test argmax(γ[:,i]) == m
        end
    end
end

@testset "GMM Classify" begin
    c1 = 6*convert(Matrix{Float64},[1 -1; 1 -1])
    c2 = 6*convert(Matrix{Float64},[1 -1; -1 1])
    μ₁ = [c1[:,i] for i ∈ axes(c1,2)]
    μ₂ = [c2[:,i] for i ∈ axes(c2,2)]
    Σ1 = [1 0; 0 1]
    Σ2 = [1 0; 0 1]
    Σ3 = [1 0; 0 1]
    Σ4 = [1 0; 0 1]
    Σ₁ = convert(Vector{Matrix{Float64}},[Σ1, Σ2])
    Σ₂ = convert(Vector{Matrix{Float64}},[Σ3, Σ4])
    w₁ = [0.5, 0.5]
    w₂ = [0.5, 0.5]

    g1 = SpeechBox.GMM(w₁,μ₁,Σ₁)
    g2 = SpeechBox.GMM(w₂,μ₂,Σ₂)

    x1 = SpeechBox.sample(g1,100000)
    x2 = SpeechBox.sample(g2,100000)

    t1 = SpeechBox.sample(g1,1000)
    t2 = SpeechBox.sample(g2,1000)
    
    gm1,ll1 = SpeechBox.trainML(SpeechBox.init_rand(), x1, 2, 10)
    gm2,ll2 = SpeechBox.trainML(SpeechBox.init_rand(), x2, 2, 10)

    for i ∈ axes(t1,2)
        @test SpeechBox.logprob(gm1,t1[:,i]) > SpeechBox.logprob(gm2,t1[:,i])
    end
    
    for i ∈ axes(t2,2)
        @test SpeechBox.logprob(gm1,t2[:,i]) < SpeechBox.logprob(gm2,t2[:,i])
    end
end