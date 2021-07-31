# Voice Activity Detection

The `vad` function takes a `framed_signal` and generates an array 1s and 0s with each element corresponding to a frame and with a 1 denoting a frame of voiced speech and a 0 denoting unvoiced speech (or an absence of speech). The `vad` function can make use different voice activity detection algorithms.

## VAD algorithms
The Voice Activity Detection algorithms currently implemented are:

- **Energy Threshold:** Assigns every frame with *energy* greater than a selected *threshold* as voiced. The *threshold* is given relative to the maximum frame energy across all input frames.
- **Energy Fraction:** assigns a certain fraction of the frames with the lowest *energy* as not voiced.

```@docs
vad
```

