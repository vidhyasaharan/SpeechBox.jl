# SpeechBox Documentation

```@docs
SpeechBox
```

The time-frequency core — the `signal`/`framed_signal`/`spectrum`/`timefreq` data
structures, signal framing, window functions, spectral analyses (`magspec`, `specgram`,
`periodogram`) and their supporting utilities — lives in
[TimeFrequencyAnalysis.jl](https://vidhyasaharan.github.io/TimeFrequencyAnalysis.jl/dev/)
and is re-exported by SpeechBox, so `using SpeechBox` provides the complete API. The name
`speech_waveform` is kept as an alias for the `signal` container. Plot recipes are
provided for signals, magnitude spectra, spectrograms and pitch overlays.

## Installation

Neither package is registered yet; add the core first, then SpeechBox:

```julia
pkg> add https://github.com/vidhyasaharan/TimeFrequencyAnalysis.jl.git
pkg> add https://github.com/unsw-edu-au/SpeechBox.jl.git
```

## Outline

SpeechBox.jl implements routines for

- [Voice Activity Detection](vad.md)
- [LPC Analyses](lpc.md)
- [Pitch Estimation](pitch.md)
- [MFCC extraction](mfcc.md)
- [K-Means Clustering](kmeans.md)
- [Utilities](utilities.md)

For framing, spectral analyses, correlations, filter responses, plotting and the data
structures, see the
[TimeFrequencyAnalysis.jl documentation](https://vidhyasaharan.github.io/TimeFrequencyAnalysis.jl/dev/).
