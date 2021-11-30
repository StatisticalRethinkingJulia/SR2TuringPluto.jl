### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 006a6698-f9dd-11ea-3f64-152c31c749ef
using Pkg, DrWatson

# ╔═╡ 0090af24-f9dd-11ea-0eb4-eb1758ba5442
begin
	using Distributions
	using StatsPlots, Plots
	using Turing
	using Logging
	using LaTeXStrings
nd

# ╔═╡ 9453eb8e-f9dc-11ea-1865-f168816cdb46
md"## Clip-02-06-07.jl"

# ╔═╡ 009dc1c8-f9dd-11ea-193f-3f1baf487786
md"## snippet 2.6"

# ╔═╡ 009f1ad2-f9dd-11ea-2abc-0579fc5057d2
@model globethrowing(w, l) = begin
    p ~ Uniform(0, 1)
    w ~ Binomial(w + l, p)
end

# ╔═╡ 00a7e61c-f9dd-11ea-028c-ddfdbdd620cb
m = globethrowing(6, 3);

# ╔═╡ 00a87104-f9dd-11ea-2ef2-bf5a7c31c607
r = quap(m)

# ╔═╡ ea69e428-f9e0-11ea-0c97-e3b70306a220
begin
	p_grid = range(0, 1, length = 20)
	prior = ones(20)
	likelihood = pdf.(Binomial.(9, p_grid), 6)
	posterior = likelihood .* prior
	posterior ./= sum(posterior)
	
	x = 0:0.01:1
	w = 6
	l = 3
	
	f = plot( x, pdf.(Normal(mean(r.coef.p), √r.vcov[1]), x ), lab="quap")
	plot!( x, pdf.(Beta( w+1 , l+1 ) , x ), lab="exact", leg=:topleft, title="n = $(w+l)")
end

# ╔═╡ 00dd7c78-f9dd-11ea-29ee-b9fb0b25bd56
md"## End clip-02-06-07t.jl"

# ╔═╡ Cell order:
# ╟─9453eb8e-f9dc-11ea-1865-f168816cdb46
# ╠═006a6698-f9dd-11ea-3f64-152c31c749ef
# ╠═0090af24-f9dd-11ea-0eb4-eb1758ba5442
# ╟─009dc1c8-f9dd-11ea-193f-3f1baf487786
# ╠═009f1ad2-f9dd-11ea-2abc-0579fc5057d2
# ╠═00a7e61c-f9dd-11ea-028c-ddfdbdd620cb
# ╠═00a87104-f9dd-11ea-2ef2-bf5a7c31c607
# ╠═ea69e428-f9e0-11ea-0c97-e3b70306a220
# ╟─00dd7c78-f9dd-11ea-29ee-b9fb0b25bd56
