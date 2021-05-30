module SpeechBox

using DSP
using FFTW
using StaticArrays
using RecipesBase
using Interpolations
using Random
using LoopVectorization
# using LinearAlgebra
# using SampledSignals

export speech_waveform, framed_signal, extract_frame, dftspec, specgram, magspec, melbankm, melfcc, frq2mel, mel2frq,
                Float, frame_energy, vad, periodogram, spectrum, timefreq, amp2db, pow2db, lpc, lpc_response,
                acf, nccf

const Float = Float64 #Set to Float32 to for 32bit floating point operations - NOT YET IMPLEMENTED or TESTED

include("signal_objects.jl")
include("internal_utilities.jl")
include("loopvectorized_utilities.jl")
include("speechframing.jl")
include("spectralanalyses.jl")
include("mfcc.jl")
include("vad.jl")
include("lpc.jl")
include("spectral_comb.jl")
include("pitch_spectral_comb.jl")
include("correlations.jl")
include("plot_recipes.jl")


end # module
