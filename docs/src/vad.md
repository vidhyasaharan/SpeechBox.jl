# `Voice Activity Detection`

The `vad` function takes a `framed_signal` and generates an array 1s and 0s with each element corresponding to a frame and with a 1 denoting a frame of voiced speech and a 0 denoting unvoiced speech (or an absence of speech). The `vad` function can make use different voice activity detection algorithms.

```@docs
vad
```

The Voice Activity Detection algorithms currently implemented are:

```@docs
SpeechBox.vad_energy_threshold
SpeechBox.vad_energy_fraction
```
