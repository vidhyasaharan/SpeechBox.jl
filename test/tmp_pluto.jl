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
	using Revise
	using SpeechBox
	using Plots; plotlyjs()
	using WAV
	using PlutoUI
end

# ╔═╡ c0757d27-7190-43ab-8273-1d9b4abbd322
md"""
# Magnitude Response
"""

# ╔═╡ 43cecbf7-d7a8-48e6-a861-5dd708cd6973
# begin
# 	p = plot(cr[:,fnum], legend=false, title="Frame Number: $fnum")
# 	plot!(p,ii,mm, seriestype=:scatter, legend = false)
# end

# ╔═╡ ab73bb62-601b-4b54-b975-f029b67dcb8b
# begin
# 	plot(cr[:,fnum], legend = false)
# 	plot!(cr[:,fnum+1], legend = false, seriescolor = :green)
# end

# ╔═╡ 887c139e-743e-4eb3-b785-7ae1d9660807
begin
	testdir = normpath(joinpath(dirname(pathof(SpeechBox)),"../test/"));
	x, fs = wavread(joinpath(testdir,"King.wav"));
	frames = framed_signal(x,fs,0.09,0.01);
end

# ╔═╡ ced76953-cb83-4fda-a9e5-14bc3867d983
@bind fnum Slider(1:frames.num_signal_frames-1)

# ╔═╡ ff88d51b-537b-4eb4-948b-297ae2008d52
delta = SpeechBox.comb_resp_Δ(frames);

# ╔═╡ a800d3a9-06fb-47a4-b33e-0493a69f6953
begin
	plot(exp.(-delta*10^7), legend = false)
	vline!([fnum])
end

# ╔═╡ 0db58d76-5e51-4717-9f4e-6c902573390b
prob = SpeechBox.prob_voiced(frames);

# ╔═╡ 75935244-9e37-4dbe-9beb-e127292401c5
pk_diff = SpeechBox.comb_resp_peak_diff(frames);

# ╔═╡ e2ec33de-b559-448b-b110-b69ffc8f3352
pk_dist = SpeechBox.dist_pitch_range(frames);

# ╔═╡ a9e00616-c2fb-4696-8c26-a30f43b8a124
cr, frqs = SpeechBox.xcorr_spectral_comb(frames);

# ╔═╡ 98ef5279-8412-443b-a051-2b9ae765d28e
ii,mm = SpeechBox.findpeaks_sorted(cr[:,fnum]);

# ╔═╡ 7cf135cf-cd13-497c-9a14-52b5e888c5ac
frame = extract_frame(frames,fnum);

# ╔═╡ 986795ba-df19-41b8-9450-3f1625fd882e
msp = SpeechBox.specgram(frames);

# ╔═╡ 1f67d634-0f0a-469e-9e42-8b75080d8e1c
begin
	p2 = plot(log(msp))
	vline!(p2,[fnum],seriescolor=:green,legend=false)
end

# ╔═╡ Cell order:
# ╟─c0757d27-7190-43ab-8273-1d9b4abbd322
# ╠═43cecbf7-d7a8-48e6-a861-5dd708cd6973
# ╠═ab73bb62-601b-4b54-b975-f029b67dcb8b
# ╟─1f67d634-0f0a-469e-9e42-8b75080d8e1c
# ╠═ced76953-cb83-4fda-a9e5-14bc3867d983
# ╠═a800d3a9-06fb-47a4-b33e-0493a69f6953
# ╠═ff88d51b-537b-4eb4-948b-297ae2008d52
# ╠═0db58d76-5e51-4717-9f4e-6c902573390b
# ╠═75935244-9e37-4dbe-9beb-e127292401c5
# ╠═e2ec33de-b559-448b-b110-b69ffc8f3352
# ╠═a9e00616-c2fb-4696-8c26-a30f43b8a124
# ╠═98ef5279-8412-443b-a051-2b9ae765d28e
# ╠═dc66e5c6-f5d9-4637-9dae-46c6a34ac79b
# ╠═887c139e-743e-4eb3-b785-7ae1d9660807
# ╠═7cf135cf-cd13-497c-9a14-52b5e888c5ac
# ╠═986795ba-df19-41b8-9450-3f1625fd882e
