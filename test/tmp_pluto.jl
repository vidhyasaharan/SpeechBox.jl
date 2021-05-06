### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ dc66e5c6-f5d9-4637-9dae-46c6a34ac79b
begin
	using SpeechBox
	using Plots; plotly()
	using WAV
	using PlutoUI
end

# ╔═╡ c0757d27-7190-43ab-8273-1d9b4abbd322
md"""
# Magnitude Response
"""

# ╔═╡ 887c139e-743e-4eb3-b785-7ae1d9660807
testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"));

# ╔═╡ 758591e2-fb0b-4ffe-87bb-7eea84ba45e0
x, fs = wavread(joinpath(testdir,"King.wav"));

# ╔═╡ 1daae57f-f439-4f76-983b-e14a236e2d98
frames = framed_signal(x,fs,0.09,0.01);

# ╔═╡ ced76953-cb83-4fda-a9e5-14bc3867d983
@bind fnum Slider(1:frames.num_signal_frames)

# ╔═╡ 7cf135cf-cd13-497c-9a14-52b5e888c5ac
frame = extract_frame(frames,fnum);

# ╔═╡ 7cce8deb-a312-4421-8f64-f339bb3552ed
mag = SpeechBox.magspec(frame,fs);

# ╔═╡ a2ac0323-f504-47d3-95ee-4285fc4e9a1b
plot(log(mag))

# ╔═╡ 62d46765-e722-434b-bfda-3957b7b2c0c7
h = SpeechBox.lpc_response(frame,fs);

# ╔═╡ 7351e405-9c5b-411f-aa70-95d6105ec74d
plot(log(h))

# ╔═╡ 986795ba-df19-41b8-9450-3f1625fd882e
msp = SpeechBox.specgram(frames);

# ╔═╡ 0837d728-4f35-4b11-85af-f21dbd3128ce
plot(log(msp))

# ╔═╡ Cell order:
# ╠═c0757d27-7190-43ab-8273-1d9b4abbd322
# ╟─ced76953-cb83-4fda-a9e5-14bc3867d983
# ╠═a2ac0323-f504-47d3-95ee-4285fc4e9a1b
# ╠═7351e405-9c5b-411f-aa70-95d6105ec74d
# ╠═0837d728-4f35-4b11-85af-f21dbd3128ce
# ╠═7cce8deb-a312-4421-8f64-f339bb3552ed
# ╠═62d46765-e722-434b-bfda-3957b7b2c0c7
# ╠═dc66e5c6-f5d9-4637-9dae-46c6a34ac79b
# ╠═887c139e-743e-4eb3-b785-7ae1d9660807
# ╠═758591e2-fb0b-4ffe-87bb-7eea84ba45e0
# ╠═1daae57f-f439-4f76-983b-e14a236e2d98
# ╠═7cf135cf-cd13-497c-9a14-52b5e888c5ac
# ╠═986795ba-df19-41b8-9450-3f1625fd882e
