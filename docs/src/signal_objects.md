# `Data Structures`

The speech processing routines in SpeechBox.jl take as input either a speech signal or frames of a speech signal. Two structs are defined for this purpose:
 - `speech_waveform` - Speech waveform in vector and sampling rate
 - `framed_signal` - Speech waveform as a vector with sampling rate and information for framing the signal

The spectral analyses routines compute the spectral components of a signal (or frame) and the following structs are defined to store this information. Plotting recipes are provided for these structs and `plot()` from Plots.jl can be called on these structs to plot the information contained in them.
 - `spectrum` - Speech waveform along with spectral components, frequency indices corresponding to spectral components and a title
 - `timefreq` - Speech waveform with framing information, time-frequency components, frequency and time indices and a title


## `speech_waveform`
Stores a time domain waveform in a vector `x` along with the sampling rate `fs`. Most relevant functions accept a signal array and sampling rate as distinct inputs but the use of `speech_waveform` objects avoids confusion when multiple signals with different sampling rates are involved or confusion with array dimensions.

```julia
struct speech_waveform
    x::Vector{Float}
    fs::Float
end
```

A straightforward constructor is used:
```@docs
speech_waveform
```

## `framed_signal`
Stores a [`speech_waveform`] (@ref) as well as all the information required to split the signal into frames. Namely, the length of the frame (`frame_length`), the number of samples between the start of consecutive frames (`frame_shift`), the total number of frames covering the signal with the smallest amount of zero padding to ensure an integer number of frames (`num_frames`), and the total number of frames spanning as much of the signal as possible with an integer number of frames without requiring any padding and ignoring any samples at the end of the signal that do not fit into a full frame (`num_signal_frames`).

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
Stores a `speech_waveform`, an array of spectral components (`components`), a vector of frequency indices corresponding to each spectral component (`frqs`), and a title string (`title`).

```julia
struct spectrum{T<:RorC}
    signal::speech_waveform
    components::Vector{T}
    frqs::Vector{Float}
    title::AbstractString
end
```

Two constructors are provided:
```@docs
spectrum
```

## `timefreq`
Stores a `speech_waveform` along with a matrix comprising of spectro-temporal compoenents (`components`) where each column corresponds to a different time and each row to a different frequency , the `framed_signal` used to obtain the spectro-temporal decomposition (if relevant), a vector of frequency indices corresponding to each row (`frqs`), either a vector of time indices, each element corresponding to each column or a vecctor of vector of time indices if time indices are not consistent across spectral components (`time`), and a title string (`title`)

```julia
struct timefreq
    signal::speech_waveform
    frames::Union{framed_signal, Nothing}
    components::RorC_Matrix
    frqs::Vector{Float}
    time::Union{Vector{Float}, Vector{Vector{Float}}}
    title::Union{AbstractString, Nothing}
end
```

Five contructors are provided:
```@docs
timefreq
```