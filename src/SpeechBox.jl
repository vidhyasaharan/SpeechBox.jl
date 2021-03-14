module SpeechBox

using DSP
using FFTW
using StaticArrays
using RecipesBase
using Interpolations
# using LinearAlgebra
# using SampledSignals

export speech_waveform, framed_signal, extract_frame, specgram, magspec, melbankm, melfcc,
                Float, frame_energy, vad, periodogram, spectrum, timefreq, amp2db, pow2db

const Float = Float64 #Set to Float32 to for 32bit floating point operations - NOT YET IMPLEMENTED or TESTED

include("signal_objects.jl")
include("internal_utilities.jl")
include("speechframing.jl")
include("spectralanalyses.jl")
include("mfcc.jl")
include("vad.jl")
include("lpc.jl")
include("plot_recipes.jl")


end # module
