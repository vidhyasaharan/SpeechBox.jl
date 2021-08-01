# Utility Functions

Some utility functions that are commonly employed in speech processing systems and algorithms are provided:
- `preemphasis` - Pre-emphasise speech signal using high pass filter: ``H(z) = 1 - 0.95z^{-1}``
- `resample` - Resample speech signal given new sampling rate
- `window` - Generate a array corresponding to a window function

## Pre-Emphasis
```@docs
preemphasis
```

## Resample
```@docs
resample
```

## Window Function
```@docs
window
```