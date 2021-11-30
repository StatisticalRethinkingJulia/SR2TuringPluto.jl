### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 849e15e4-f786-11ea-3e89-315057d10f40
using Pkg, DrWatson

# ╔═╡ 849e4d66-f786-11ea-1fff-71edc16a1bdb
begin
	@quickactivate "SR2TuringPluto"
	using StatisticalRethinking
end

# ╔═╡ 2b094a8a-f786-11ea-3688-3769dfa0c5f8
md"## Fig3.2t.jl"

# ╔═╡ 849eb8c8-f786-11ea-1b05-e7a8c5c1d987
begin
	p_grid = range(0, stop=1, length=10000)
	prior = ones(length(p_grid))
	likelihood = pdf.(Binomial.(9, p_grid), 6)
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
	samples = sample(p_grid, Weights(posterior), length(p_grid));
end;

# ╔═╡ 84abc770-f786-11ea-3a46-b5f847b0f5b8
begin
	N = 10000
	samples2 = sample(p_grid, Weights(posterior), N)
end;

# ╔═╡ 84ac5e4c-f786-11ea-37bc-3318feecfff9
begin
	b1 = mapreduce(p -> p < 0.5 ? 1 : 0, +, samples2) / N
	b2 = mapreduce(p -> (p > 0.5 && p < 0.75) ? 1 : 0, +, samples2) / N
	b3 = quantile(samples2, 0.8)
	b4 = quantile(samples2, [0.1, 0.9])

	p1 = plot_density_interval(samples2, [0.0, 0.5],
	  xlab="Proportion water (p)");
	p2 = plot_density_interval(samples2, [0.5, 0.75],
	  xlab="Proportion water (p)");
	p3 = plot_density_interval(samples2, [0.0, b3], 
	  xlab="Proportion water (p)");
	p4 = plot_density_interval(samples2, b4, 
	  xlab="Proportion water (p)")
end;

# ╔═╡ 84ba27de-f786-11ea-00e5-ff77781b2b88
plot(p1, p2, p3, p4, layout=(2, 2))

# ╔═╡ 84c23eb0-f786-11ea-262c-2f906a343b13
md"## End of Fig3.2t.jl"

# ╔═╡ Cell order:
# ╟─2b094a8a-f786-11ea-3688-3769dfa0c5f8
# ╠═849e15e4-f786-11ea-3e89-315057d10f40
# ╠═849e4d66-f786-11ea-1fff-71edc16a1bdb
# ╠═849eb8c8-f786-11ea-1b05-e7a8c5c1d987
# ╠═84abc770-f786-11ea-3a46-b5f847b0f5b8
# ╠═84ac5e4c-f786-11ea-37bc-3318feecfff9
# ╠═84ba27de-f786-11ea-00e5-ff77781b2b88
# ╟─84c23eb0-f786-11ea-262c-2f906a343b13
