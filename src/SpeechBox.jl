"""
    SpeechBox

Basic speech processing and analysis routines: framing-based feature extraction (MFCCs),
voice activity detection, LPC analyses, pitch estimation, and clustering / Gaussian mixture
modelling utilities.

The time-frequency core (signal containers, framing, windows, spectral analyses,
correlations, filter responses and their supporting utilities) lives in
[TimeFrequencyAnalysis.jl](https://github.com/vidhyasaharan/TimeFrequencyAnalysis.jl) and is
re-exported here, so `using SpeechBox` provides the complete API.
"""
module SpeechBox

using Reexport
@reexport using TimeFrequencyAnalysis #time-frequency core: signal/framed_signal/spectrum/timefreq, framing, windows, spectral analyses, correlations, filter responses, grids, generators
using TimeFrequencyAnalysis: findclosest, generate_ticks #unexported helpers: k-means++/categorical sampling; axis ticks for the pitch plot recipes

import DSP #qualified access to DSP.Filters (vocal tract models in lpc.jl)
using DSP: hamming, nextfastfft
using FFTW: plan_rfft, plan_dct
using LinearAlgebra
using RecipesBase

#The signal container was named speech_waveform before the core moved to
#TimeFrequencyAnalysis; the historical name is kept as an alias for backwards
#compatibility
const speech_waveform = TimeFrequencyAnalysis.signal
export speech_waveform

export preemphasis

export vad, energy_threshold, energy_fraction

export lpc, lpc_response

export melbankm, melfcc, frq2mel, mel2frq

export pitch, spectral_comb

export kmeans, kmeans!, kmpp, kmrand, kmeans_init, closest_centre, closest_centre!, mindist2cntrs

include("pitch_objects.jl")
include("utilities.jl")
include("internal_utilities.jl")
include("statistics_utilities.jl")
include("mfcc.jl")
include("vad.jl")
include("lpc.jl")
include("spectral_comb.jl")
include("pitch.jl")
include("pitch_RAPT.jl")
include("plot_recipes.jl")
include("levinson_durbin.jl")
include("kmeans.jl")
include("GMM.jl")

end # module
