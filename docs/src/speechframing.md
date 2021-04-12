# `Speech Frames`

Two basic functions are provided to manipulate speech frames:
- `extract_frame` - Extract desired frame from a [`framed_signal`](@ref) object
- `frame_energy` - Estimate signal 'energy' in each frame

```@docs
lpc
```

## `extract_frame`
The primary specech framing function is [`extract_frame`](@ref) which is used to extract the desired frame from a speech signal.

```@docs
extract_frame
```

## `frame_energy`

The utility function [`frame_energy`](@ref) computes the 'energy' (L2 norm) in each frame of a signal. Optionally, the output can be 'normalised' such that the maximum frame energy is 1.

```@docs
frame_energy
```