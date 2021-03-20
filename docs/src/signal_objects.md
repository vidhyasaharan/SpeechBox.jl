# `Signal Objects`

SpeechBox.jl defines 4 structs:
 - `speech_waveform` - Speech waveform in vector and sampling rate
 - `framed_signal` - Speech waveform as a vector with sampling rate and information for framing the signal
 - `spectrum` - Speech waveform along with spectral components, frequency indices corresponding to spectral components and a title
 - `timefreq` - Speech waveform with framing information, time-frequency components, frequency and time indices and a title


## speech_waveform
The struct stores a time domain waveform in a vector `x` along with the sampling rate `fs`. Most relevant functions accept a signal array and sampling rate as distinct inputs but the use of `speech_waveform` object avoids confusion when multiple signals with different sampling rates are involved or confusion with array dimensions.

```julia
struct speech_waveform
    x::Array{Float,1}
    fs::Float
end
```



```@docs
speech_waveform
framed_signal
spectrum
timefreq
```