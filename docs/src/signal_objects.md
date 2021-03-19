# `Signal Objects`

SpeechBox.jl defines 4 datatypes
 - `speech_waveform` - Speech waveform in vector and sampling rate
 - `framed_signal` - Speech waveform as a vector with sampling rate and information for framing the signal
 - `spectrum` - Speech waveform along with spectral components, frequency indices corresponding to spectral components and a title
 - `timefreq` - Speech waveform with framing information, time-frequency components, frequency and time indices and a title

```@docs
speech_waveform
framed_signal
spectrum
timefreq
```