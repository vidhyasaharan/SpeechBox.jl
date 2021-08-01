# Pitch Estimation
As is common usage, the term pitch estimation actually refers to estimation of the (time-varying) fundamental frequency. The following methods/algorithms for pitch estimation are implemented:
- *Spectral Comb* - A VAD is used to identy voiced frames and the pitch (fundamental frequency) in each voiced frame is estimated based on the position of the peak of the response to a spectral comb.
- *RAPT* - Not implemented yet

```@docs
pitch
```
