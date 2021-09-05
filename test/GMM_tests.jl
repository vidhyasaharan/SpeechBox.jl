data = SpeechBox.generate_4_clusters(5)
g = SpeechBox.GMMinit(4,data)

p = scatter(data[1,:], data[2,:])
scatter!(p,g.μ[1])
scatter!(p,g.μ[2])
scatter!(p,g.μ[3])
scatter!(p,g.μ[4])



x = randn(10,100000)
m1 = Vector{Float64}(undef,size(x,1))
v1 = Vector{Float64}(undef,size(x,1))
m = Vector{Float64}(undef,size(x,1))
C = Matrix{Float64}(undef,size(x,1),size(x,1))
@benchmark SpeechBox.running_meanvar!(m1,v1,x)
@benchmark SpeechBox.running_meancov!(m,C,x)