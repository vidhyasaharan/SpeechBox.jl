# `Signal Objects`

SpeechBox.jl defines 4 structs:
 - `speech_waveform` - Speech waveform in vector and sampling rate
 - `framed_signal` - Speech waveform as a vector with sampling rate and information for framing the signal
 - `spectrum` - Speech waveform along with spectral components, frequency indices corresponding to spectral components and a title
 - `timefreq` - Speech waveform with framing information, time-frequency components, frequency and time indices and a title


## `speech_waveform`
Stores a time domain waveform in a vector `x` along with the sampling rate `fs`. Most relevant functions accept a signal array and sampling rate as distinct inputs but the use of `speech_waveform` objects avoids confusion when multiple signals with different sampling rates are involved or confusion with array dimensions.

```julia
struct speech_waveform
    x::Array{Float,1}
    fs::Float
end
```

A straightforward constructor is used:
```@docs
speech_waveform
```

## `framed_signal`
Stores a [`speech_waveform`] (@ref `speech_waveform`) as well as all the information required to split the signal into frames. Namely, the length of the frame (`frame_length`), the number of samples between the start of consecutive frames (`frame_shift`), the total number of frames covering the signal with the smallest amount of zero padding to ensure an integer number of frames (`num_frames`), and the total number of frames spanning as much of the signal as possible with an integer number of frames without requiring any padding and ignoring any samples at the end of the signal that do not fit into a full frame (`num_signal_frames`).

```julia
struct framed_signal
    signal::speech_waveform
    frame_length::Int
    frame_shift::Int
    num_signal_frames::Int
    num_frames::Int
end
```

Two constructors are provided:
```@docs
framed_signal
```

## `spectrum`
Stores