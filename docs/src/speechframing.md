# Speech Frames

Four basic functions are provided to manipulate speech frames:
- `view_frame` - Create a [view](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-views) of the desired frame from a [`framed_signal`](@ref) object
- `extract_frame` - Extract desired frame from a [`framed_signal`](@ref) object
- `enframe` - Store all frames of a signal as columns in a matrix.
- `frame_energy` - Estimate signal 'energy' in each frame

## `view_frame`
The primary speech framing function is [`view_frame`](@ref) which is used to create a [view](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-views) of the desired frame from a speech signal.

```@docs
view_frame
```

## `extract_frame`
The speech framing function [`extract_frame`](@ref) is similar to `view_frame`, but is extracts the desired frame from a speech signal as a vector instead of creating a view.

```@docs
extract_frame
```

## `enframe`
The function [`enframe`](@ref) stores frames of a signal as columns of a matrix.

```@docs
enframe
enframe!
```


## `frame_energy`

The utility function [`frame_energy`](@ref) computes the 'energy' (L2 norm) in each frame of a signal. Optionally, the output can be 'normalised' such that the maximum frame energy is 1.

```@docs
frame_energy
```