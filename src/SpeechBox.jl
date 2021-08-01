module SpeechBox


using DSP
using FFTW
using RecipesBase
using Random
using LoopVectorization



export speech_waveform, framed_signal, spectrum, timefreq

export preemphasis, window, resample

export extract_frame, view_frame, enframe, enframe!, frame_energy

export dft, magspec, specgram, periodogram

export vad, energy_threshold, energy_fraction

export lpc, lpc_response

export xcorr!, xcorr, acorr!, acorr, acf, nccf

export melbankm, melfcc, frq2mel, mel2frq

export pitch, spectral_comb

export amp2db, pow2db

export comp

const Float = Float64 #Set to Float32 to for 32bit floating point operations - NOT YET IMPLEMENTED or TESTED

struct comp end

include("signal_objects.jl")
include("utilities.jl")
include("internal_utilities.jl")
include("loopvectorized_utilities.jl")
include("dsp_utilities.jl")
include("speechframing.jl")
include("spectralanalyses.jl")
include("mfcc.jl")
include("vad.jl")
include("lpc.jl")
include("spectral_comb.jl")
include("pitch.jl")
include("correlations.jl")
include("pitch_RAPT.jl")
include("plot_recipes.jl")
include("levinson_durbin.jl")
include("distortion_measures.jl")


end # module
