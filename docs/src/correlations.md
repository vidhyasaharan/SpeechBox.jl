# Correlation

Functions for autocorrelation and cross-correlation of discrete-time sequences, typically frames of speech, are provided. Additionally, autocorrelation function (ACF) and normalised autocorrelation function (NACF) implementation intended for entire speech segments are also provided.

## Cross-Correlation Sequences
The [`xcorr`](@ref) and [`xcorr!`](@ref) functions compute the cross correlation between two sequences and would typically be used to search along a longer sequence for the best match to the shorter sequence. It is implemented as a *sliding dot product*.

```@docs
xcorr
xcorr!
```

## Autocorrelation Sequence
The [`acorr`](@ref) and [`acorr!`](@ref) functions compute autocorrelation sequence up to the desired lag. The autocorrelation value at each lag is normalised by the size of the correlation window, which becomes smaller at the edges.

```@docs
acorr
acorr!
```

## Autocorrelation Functions
The [`acf`](@ref) function returns autocorrelation function value of a signal at a specified position, at one or more lags, using a specified window size and type. The normalised version, [`nacf`](@ref), returns the autocorrelation computed using windows normalised by mean subtraction and scaled by their norm.

```@docs
acf
nacf
```
