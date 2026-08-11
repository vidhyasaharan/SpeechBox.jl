# LPC Analyses

Linear Predictive Coding (LPC) is based on the linear all-pole vocal tract filter model. It assumes short frames of speech can be modelled as an autoregressive (AR) process. The coefficients of the corresponding all-pole model (AR model parameters) are inferred using the Levinson-Durbin method. The following functions are provided:

- `lpc` - Estimate the AR filter model coefficients given a signal
- `lpc_response` - Compute the magnitude respone of the all-pole AR filter that fits a given signal

```@docs
lpc
lpc_response
```

The autocorrelation (`acorr`) and filter-response (`filter_coefs`, `filter_magresp`) routines used in LPC analysis are provided and documented by [TimeFrequencyAnalysis.jl](https://vidhyasaharan.github.io/TimeFrequencyAnalysis.jl/dev/).