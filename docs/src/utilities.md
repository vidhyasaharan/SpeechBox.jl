# Utility Functions

Speech-specific signal conditioning provided by SpeechBox:

- `preemphasis` - Pre-emphasise speech signal using high pass filter: ``H(z) = 1 - 0.95z^{-1}``

Generic signal utilities — window functions, resampling, frequency grids and synthetic
signal generators — are provided (and documented) by
[TimeFrequencyAnalysis.jl](https://unsw-edu-au.github.io/TimeFrequencyAnalysis.jl/dev/)
and re-exported by SpeechBox.

## Pre-Emphasis
```@docs
preemphasis
```
