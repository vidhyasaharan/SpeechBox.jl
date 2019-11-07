module SpeechBox


using DSP
using FFTW
using StaticArrays
# using SampledSignals

export framed_signal, extract_frame, specgram, magspec, melbankm, melfcc, Float, frame_energy, vad

const Float = Float64 #Set to Float32 to fore 32bit floating point operations - NOT YET IMPLEMENTED or TESTED

include("internal_utilities.jl")
include("speechframing.jl")
include("spectralanalyses.jl")
include("mfcc.jl")
include("vad.jl")


end # module
