# SpeechBox

<!-- [![codecov](https://codecov.io/gh/vidhyasaharan/SpeechBox.jl/branch/master/graph/badge.svg?token=SOB3LPJO8I)](https://codecov.io/gh/vidhyasaharan/SpeechBox.jl) -->

[![CI](https://github.com/unsw-edu-au/SpeechBox.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/unsw-edu-au/SpeechBox.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://unsw-edu-au.github.io/SpeechBox.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://unsw-edu-au.github.io/SpeechBox.jl/dev)


SpeechBox is a speech processing toolbox consisting of Julia routines that are maintained by and mostly written by Dr. Vidhyasaharan Sethu, Speech and Behavioural Signal Processing Lab, School of EE&T, UNSW Australia.

SpeechBox builds on [TimeFrequencyAnalysis.jl](https://github.com/unsw-edu-au/TimeFrequencyAnalysis.jl), which provides the signal containers, framing and spectral analyses; SpeechBox re-exports that API, so `using SpeechBox` provides everything.

## Installation

Neither package is registered yet; add the core first, then SpeechBox:

```julia
pkg> add https://github.com/unsw-edu-au/TimeFrequencyAnalysis.jl.git
pkg> add https://github.com/unsw-edu-au/SpeechBox.jl.git
```
