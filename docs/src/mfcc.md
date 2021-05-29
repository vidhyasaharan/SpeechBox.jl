# Mel Frequency Cepstral Coefficients

Mel Frequency Cepstral Coefficients (MFCCs) are extracted on a frame by frame basis, with one set of coefficients per frame. The extraction is implemented in the frequency domain using the DFT implementation from [FFTW.jl](https://github.com/JuliaMath/FFTW.jl). MFCC extraction is provided as one function:
- `melfcc`- Extract a MFCC vector per frame of the input speech signal

```@docs
melfcc
melbankm
frq2mel
mel2frq
```

