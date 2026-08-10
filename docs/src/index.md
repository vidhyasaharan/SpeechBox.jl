# SpeechBox Documentation

```@docs
SpeechBox
```

The time-frequency core — the `waveform`/`framed_signal`/`spectrum`/`timefreq` data
structures, signal framing, window functions, spectral analyses (`magspec`, `specgram`,
`periodogram`) and their supporting utilities — lives in
[TimeFrequencyAnalysis.jl](https://unsw-edu-au.github.io/TimeFrequencyAnalysis.jl/dev/)
and is re-exported by SpeechBox, so `using SpeechBox` provides the complete API. The name
`speech_waveform` is kept as an alias for the `waveform` container. Plot recipes are
provided for waveforms, magnitude spectra, spectrograms and pitch overlays.

## Installation

Neither package is registered yet; add the core first, then SpeechBox:

```julia
pkg> add https://github.com/unsw-edu-au/TimeFrequencyAnalysis.jl.git
pkg> add https://github.com/unsw-edu-au/SpeechBox.jl.git
```

## Outline

SpeechBox.jl implements routines for

- [Voice Activity Detection](vad.md)
- [Correlations](correlations.md)
- [LPC Analyses](lpc.md)
- [Pitch Estimation](pitch.md)
- [MFCC extraction](mfcc.md)
- [K-Means Clustering](kmeans.md)
- [Utilities](utilities.md)

For framing, spectral analyses and the data structures, see the
[TimeFrequencyAnalysis.jl documentation](https://unsw-edu-au.github.io/TimeFrequencyAnalysis.jl/dev/).
