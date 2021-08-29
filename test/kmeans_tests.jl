ndim = 20
npts = 10000
x = randn(Float64,ndim,npts)



a = randn(Float64,10)
b = randn(Float64,10)
c = [a b a b]
@benchmark d1 = sum(abs2,a-b)
@benchmark d2 = SpeechBox.dotavx(a-b)
@benchmark d3 = SpeechBox.sqL2avx(a,b)
@benchmark d4 = SpeechBox.L2avx(a,b)
@benchmark dd = SpeechBox.pairwise(SpeechBox.SqL2(), a, c)
@benchmark dd2 = SpeechBox.pairwise(SpeechBox.SqL2(), c, c)