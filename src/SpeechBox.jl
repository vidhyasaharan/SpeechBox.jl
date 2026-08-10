module SpeechBox


using DSP
using FFTW
using RecipesBase
using Random
using LinearAlgebra



export speech_waveform, framed_signal, spectrum, timefreq

export preemphasis, window, resample

export extract_frame, view_frame, enframe, enframe!, frame_energy

export dft, magspec, specgram, periodogram

export vad, energy_threshold, energy_fraction

export lpc, lpc_response

export xcorr!, xcorr, acorr!, acorr, acf, nacf

export melbankm, melfcc, frq2mel, mel2frq

export pitch, spectral_comb

export amp2db, pow2db

export comp

export kmeans, kmeans!, kmpp, kmrand, kmeans_init, closest_centre, closest_centre!, mindist2cntrs

struct comp end

include("signal_objects.jl")
include("utilities.jl")
include("internal_utilities.jl")
include("statistics_utilities.jl")
include("dsp_utilities.jl")
include("speechframing.jl")
include("spectralanalyses.jl")
include("mfcc.jl")
include("vad.jl")
include("lpc.jl")
include("spectral_comb.jl")
include("pitch.jl")
include("pitch_RAPT.jl")
include("correlations.jl")
include("plot_recipes.jl")
include("levinson_durbin.jl")
include("kmeans.jl")
include("GMM.jl")



end # module
