### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 9e38d21a-f785-11ea-1ce0-5d313a4e9a7e
using Pkg, DrWatson

# ╔═╡ 9e393624-f785-11ea-077b-cb4f64059bb3
begin
	@quickactivate "SR2TuringPluto"
	using StatisticalRethinking
end

# ╔═╡ 4e27a698-f785-11ea-1ded-a548d1953a14
md"## Fig3.1t.jl"

# ╔═╡ 9e3b92a2-f785-11ea-090f-a78cbded0a8d
begin
	N = 10000
	p_grid = range(0, stop=1, length=N)
	prior = ones(length(p_grid))
	likelihood = pdf.(Binomial.(9, p_grid), 6)
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
	samples = sample(p_grid, Weights(posterior), length(p_grid));
end;

# ╔═╡ 9e43adc0-f785-11ea-2b6a-f11004009094
begin
	p1 = scatter(samples, ylim=(0, 1), xlab="Sample number",
	  ylab="Proportion water(p)", leg=false)
	p2 = density(samples, xlim=(0.0, 1.0), ylim=(0.0, 3.0),
	  xlab="Proportion water (p)",
	  ylab="Density", leg=false)
end;

# ╔═╡ d13a98ea-f785-11ea-0cb5-e51b20ddcdac
plot(p1, p2, layout=(1,2))

# ╔═╡ 9e520366-f785-11ea-1d50-8f7ab0c29c9c
md"## End of Fig3.1t.jl"

# ╔═╡ Cell order:
# ╟─4e27a698-f785-11ea-1ded-a548d1953a14
# ╠═9e38d21a-f785-11ea-1ce0-5d313a4e9a7e
# ╠═9e393624-f785-11ea-077b-cb4f64059bb3
# ╠═9e3b92a2-f785-11ea-090f-a78cbded0a8d
# ╠═9e43adc0-f785-11ea-2b6a-f11004009094
# ╠═d13a98ea-f785-11ea-0cb5-e51b20ddcdac
# ╟─9e520366-f785-11ea-1d50-8f7ab0c29c9c
