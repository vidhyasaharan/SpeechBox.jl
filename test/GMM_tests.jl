ndim = 100
a = randn(Float64,ndim,ndim)
A = a'*a + LinearAlgebra.I(ndim)
x = randn(ndim)
y = randn(ndim)


x = randn(3,10)
m = Vector{Float64}(undef,size(x,1))
C = Matrix{Float64}(undef,size(x,1),size(x,1))
SpeechBox.running_meancov!(m, C, x)



x = randn(10,100000)
m1 = Vector{Float64}(undef,size(x,1))
v1 = Vector{Float64}(undef,size(x,1))
m = Vector{Float64}(undef,size(x,1))
C = Matrix{Float64}(undef,size(x,1),size(x,1))
@benchmark SpeechBox.running_meanvar!(m1,v1,x)
@benchmark SpeechBox.running_meancov!(m,C,x)