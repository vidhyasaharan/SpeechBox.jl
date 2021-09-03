using SpeechBox
using BenchmarkTools
using Plots; plotlyjs()

x = SpeechBox.generate_4_clusters(10)
# icn = SpeechBox.kmeans_init(SpeechBox.kmpp(),x,4)
# cn = SpeechBox.kmeans_init(SpeechBox.kmrand(),x,4)

# cn = copy(icn)
cn = SpeechBox.kmeans(x,4)
# @benchmark SpeechBox.kmeans!(cn,x)
