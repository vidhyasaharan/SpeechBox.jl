# `Speech Frames`


## `extract_frame`
The primary specech framing function is [`extract_frame`](@ref) which is used to extract the desired frame from a speech signal.

```@docs
extract_frame
```

## `frame_energy`

The utility function [`frame_energy`](@ref) computes the 'energy' (L2 norm) in each frame of a signal. Optionally, the output can be 'normalised' such that the maximum frame energy is 1.