### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ bab5e792-0708-11eb-2adb-b50a08f43768
using Pkg, DrWatson

# ╔═╡ bab64fe8-0708-11eb-18de-8f1741de4bda
begin
	@quickactivate "SR2TuringPluto"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ b6c9d85a-0708-11eb-0ea3-e7f80bf70fce
md"## Clip-04-37-41t.jl"

# ╔═╡ bab70528-0708-11eb-246f-15df8b5370f7
#include(srcdir("tools.jl"))

md"### snippet 4.26"

# ╔═╡ bacb3f84-0708-11eb-375c-974f4c8e6453
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)
end;

# ╔═╡ bacbe5a6-0708-11eb-15d0-bb7b033210af
md"### snippet 4.37"

# ╔═╡ badb3b82-0708-11eb-0d65-23f3b0ca7256
scatter(df.weight, df.height; lab="Observations", leg=:bottomright)

# ╔═╡ bade45fe-0708-11eb-1ab9-2f0b4f7e7331
md"### snippet 4.38"

# ╔═╡ baef1026-0708-11eb-36a2-e58bf8132c67
begin
	N = 100
	a1 = rand(Normal(178.0, 20.0), N)
	b1 = rand(Normal(0.0, 10.0), N)
end;

# ╔═╡ baf373a0-0708-11eb-3b1d-c72130f8a095
md"### snippet 4.39"

# ╔═╡ bafe27a0-0708-11eb-18a2-bb1f5950bf91
begin
	plot(xlims = extrema(df.weight), ylims = (-100, 400), xlable = "weight", ylabel = "height")
	hline!([0.0, 272.0])
	for (a, b) in zip(a1, b1)
		plot!(x, a .+ b .* (x .- x̄), color = "black", alpha = 0.2)
	end
	plot!(legend = false)
end

# ╔═╡ bb109a16-0708-11eb-1154-b19c10620afe
md"### snippet 4.40"

# ╔═╡ bb139e1e-0708-11eb-392f-278ee5791f8d
begin
	b = rand(LogNormal(0.0, 1.0), 10_000)
	density(b, xlims = (0, 5))
end

# ╔═╡ bb1dc006-0708-11eb-01bf-99a9a54e165e
md"### snippet 4.41"

# ╔═╡ bb28a7a0-0708-11eb-2b42-e185aea41459
begin
	a2 = rand(Normal(178.0, 20.0), N)
	b2 = rand(LogNormal(0.0, 1.0), N)
	plot(xlims = extrema(df.weight), ylims = (-100, 400), xlable = "weight", ylabel = "height")
	hline!([0.0, 272.0])
	foreach(zip(a2, b2)) do (a, b)
		plot!(x, a .+ b .* (x .- x̄), color = "black", alpha = 0.2)
	end
	plot!(legend = false)
end

# ╔═╡ 056edb20-070b-11eb-0150-416bc812a4ca
md"## End of clip-04-37-41t.jl"

# ╔═╡ Cell order:
# ╠═b6c9d85a-0708-11eb-0ea3-e7f80bf70fce
# ╠═bab5e792-0708-11eb-2adb-b50a08f43768
# ╠═bab64fe8-0708-11eb-18de-8f1741de4bda
# ╟─bab70528-0708-11eb-246f-15df8b5370f7
# ╠═bacb3f84-0708-11eb-375c-974f4c8e6453
# ╟─bacbe5a6-0708-11eb-15d0-bb7b033210af
# ╠═badb3b82-0708-11eb-0d65-23f3b0ca7256
# ╟─bade45fe-0708-11eb-1ab9-2f0b4f7e7331
# ╠═baef1026-0708-11eb-36a2-e58bf8132c67
# ╟─baf373a0-0708-11eb-3b1d-c72130f8a095
# ╠═bafe27a0-0708-11eb-18a2-bb1f5950bf91
# ╠═bb109a16-0708-11eb-1154-b19c10620afe
# ╠═bb139e1e-0708-11eb-392f-278ee5791f8d
# ╠═bb1dc006-0708-11eb-01bf-99a9a54e165e
# ╠═bb28a7a0-0708-11eb-2b42-e185aea41459
# ╟─056edb20-070b-11eb-0150-416bc812a4ca
