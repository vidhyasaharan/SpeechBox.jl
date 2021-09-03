x = SpeechBox.generate_4_clusters(10)
cin = SpeechBox.kmeanspp(x,4)

cc1 = Vector{Int}(undef,size(x,2))
cc = SpeechBox.closest_centre(x[:,cin],x)
SpeechBox.closest_centre!(cc1,x[:,cin],x)
