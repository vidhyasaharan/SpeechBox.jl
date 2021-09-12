@testset "normaliseWeights" begin
    nmix = 10
    w = rand(Float,nmix)
    ŵ = copy(w)
    SpeechBox.normaliseWeights!(w)
    @test sum(w) ≈ one(Float)
    for i ∈ eachindex(w,ŵ)
        @test w[i]/sum(w) ≈ ŵ[i]/sum(ŵ)
    end
end

@testset "isconsistentX" begin
    function posdefmatrix(ndim::Int)
        a = randn(Float,ndim,ndim)
        A = a*a' + LinearAlgebra.I(ndim)
        return A
    end
    nmix = 5
    ndim = 3
    w = rand(nmix)
    SpeechBox.normaliseWeights!(w)
    μ = [randn(Float,ndim) for i ∈ 1:nmix]
    Σ = [posdefmatrix(ndim) for i ∈ 1:nmix]

    @test SpeechBox.isconsistentComponents(w,μ,Σ) == true
    @test SpeechBox.isconsistentComponents(w[1:nmix-1],μ,Σ) == false
    @test SpeechBox.isconsistentComponents(w,μ[1:nmix-1],Σ) == false
    @test SpeechBox.isconsistentComponents(w,μ,Σ[1:nmix-1]) == false

    @test SpeechBox.isconsistentWeights(w) == true
    @test SpeechBox.isconsistentWeights(ones(nmix)) == false

    @test SpeechBox.isconsistentDimensions(μ,Σ) == true
    μ̄ = copy(μ)
    μ̄[nmix] = randn(Float,ndim-1)
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
    function posdefmatrix(ndim::Int)
        a = randn(Float,ndim,ndim)
        A = a*a' + LinearAlgebra.I(ndim)
        return A
    end
    ndim = 10
    nmix = 5
    w = rand(nmix)
    SpeechBox.normaliseWeights!(w)
    μ = [randn(Float,ndim) for i ∈ 1:nmix]
    Σ = [posdefmatrix(ndim) for i ∈ 1:nmix]

    G = SpeechBox.GMM(w,μ,Σ)

    @test typeof(G.w) == Vector{Float}
    @test typeof(G.μ) == Vector{Vector{Float}}
    @test typeof(G.Σ) == Vector{Matrix{Float}}
    @test typeof(G.P) == Vector{Matrix{Float}}
    @test typeof(G.Z) == Vector{Float}

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

    @test typeof(G.w) == Vector{Float}
    @test typeof(G.μ) == Vector{Vector{Float}}
    @test typeof(G.Σ) == Vector{Matrix{Float}}
    @test typeof(G.P) == Vector{Matrix{Float}}
    @test typeof(G.Z) == Vector{Float}

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
    G = SpeechBox.GMMinit(4,x)
    @test G.w == (1/4)*ones(4)
    for i = 1:4
        @test SpeechBox.mindist2cntrs(kmeans,hcat(G.μ...)) == zeros(4)
    end
end

@testset "logsumexp" begin
    nprbs = 10
    probs = rand(Float,nprbs)
    lprobs = log.(probs)
    @test SpeechBox.logsumexp(lprobs) ≈ log(sum(probs))
end

@testset "Categorical" begin
    ncat = 10
    p = rand(Float,ncat)
    SpeechBox.normaliseWeights!(p)
    c = SpeechBox.pdist2cdist(p)
    for i ∈ eachindex(p,c)
        @test c[i] ≈ sum(p[1:i])
    end
end