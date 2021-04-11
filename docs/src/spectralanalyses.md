# `Spectral Analyses`

Both spectral and Spectro-temporal analyses routines are included and listed below based on analyses techniques:
- *Fourier Spectrum* - Spectral analyses based on DFT (Discrete Fourier Transform)
- *STFT Spectrum* - Spectro-temporal analyses based on STFT (Short Time Fourier Transform)
- *Periodograms* - Spectral and Spectro-Temporal analyses based on periodograms (inner products with complex exponentials)
- *LPC Models* - Spectral and Spectro-Temporal analyses based on Linear Predictive Coding (LPC)/Autoregressive (AR) filter models of speech

## Fourier Spectrum (Spectral)
```@docs
dftspec
magspec
```

## STFT Spectrogram (Spectro-Temporal)spectra analyses docs
```@docs
specgram
```

## Periodogram (Spectral and Spectro-Temporal)
```@docs
periodogram
```

## LPC/AR Filter Model Response (Spectral and Spectro-temporal)
```@docs
lpc_response
```