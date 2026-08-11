@testset "mindist2cntrs" begin
    npts_per_cluster = 4
    r = 0.25
    c = convert(Matrix{Float64},[1 1 -1 -1; 1 -1 1 -1])
    x = SpeechBox.generate_4_circle_clusters(npts_per_cluster; r)
    md = SpeechBox.mindist2cntrs(c,x)
    for i ∈ eachindex(md)
        @test md[i] == r^2
    end
end


@testset "closest_centre" begin
    npts_per_cluster = 4
    r = 0.25
    c = convert(Matrix{Float64},[1 1 -1 -1; 1 -1 1 -1])
    x = SpeechBox.generate_4_circle_clusters(npts_per_cluster; r)
    cc = SpeechBox.closest_centre(c,x)
    for i ∈ axes(c,2)
        sindx = (i-1)*npts_per_cluster+1
        eindx = i*npts_per_cluster
        @test cc[sindx:eindx]==ones(npts_per_cluster)*i
    end
end


@testset "max_centre_shift" begin
    c = convert(Matrix{Float64},[1 1 -1 -1; 1 -1 1 -1])
    ncntrs = size(c,2)
    r = rand(ncntrs)
    θ = rand(ncntrs)*2π
    shift = Matrix{Float64}(undef,2,ncntrs)
    for i ∈ axes(c,2)
        shift[1,i] = r[i]*cos(θ[i])
        shift[2,i] = r[i]*sin(θ[i])
    end
    ĉ = c .+ shift
    mshift = SpeechBox.max_centre_shift(ĉ,c)
    @test mshift ≈ maximum(r)
end


@testset "kmeans_init" begin
    npts_per_cluster = 4
    r = 0.25
    x = SpeechBox.generate_4_circle_clusters(npts_per_cluster; r)
    
    cntrs = SpeechBox.kmeans_init(SpeechBox.kmpp(), x, 4)
    @test size(cntrs,1) == size(x,1)
    @test size(cntrs,2) == 4
    @test eltype(cntrs) == eltype(x)
    md = SpeechBox.mindist2cntrs(x,cntrs)
    @test md == zeros(4)

    cntrs = SpeechBox.kmeans_init(SpeechBox.kmrand(), x, 4)
    @test size(cntrs,1) == size(x,1)
    @test size(cntrs,2) == 4
    @test eltype(cntrs) == eltype(x)
    md = SpeechBox.mindist2cntrs(x,cntrs)
    @test md == zeros(4)
end


@testset "k-means" begin
    npts_per_cluster = 8
    r = 0.25
    x = SpeechBox.generate_4_circle_clusters(npts_per_cluster; r)
    cn = SpeechBox.kmeans_init(SpeechBox.kmpp(), x, 4)
    icn = copy(cn)
    cc = Vector{Int}(undef,size(x,2))
    SpeechBox.kmeans_update!(cc,cn,x)
    for i ∈ axes(cn,2)
        @test SpeechBox.mindist2cntrs(icn,cn[:,[i]]) > [0.0]
    end

    kcn = SpeechBox.kmeans(x,4)
    @test size(kcn,1) == size(x,1)
    @test size(kcn,2) == 4
    @test eltype(kcn) == eltype(x)
    SpeechBox.closest_centre!(cc,kcn,x)
    for i ∈ axes(kcn,2)
        @test kcn[:,i] == SpeechBox.running_mean(x[:,cc.==i])
    end
end