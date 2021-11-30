### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 15b6d7c8-f787-11ea-0fd9-4738ba47533d
using Pkg, DrWatson

# ╔═╡ 15b712ce-f787-11ea-245f-55620c5d1b05
begin
	@quickactivate "SR2TuringPluto"
	using StatisticalRethinking
end

# ╔═╡ bd4d25cc-f786-11ea-18fd-3326d3896012
md"## Fig3.3t.jl"

# ╔═╡ 15b79bfe-f787-11ea-010f-5d6a113c0f19
begin
	p_grid = range(0, step=0.000001, stop=1)
	prior = ones(length(p_grid))
	likelihood = pdf.(Binomial.(3, p_grid), 3)
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
	samples = sample(p_grid, Weights(posterior), length(p_grid));
end;

# ╔═╡ 15bf2d08-f787-11ea-1316-e71570fded77
begin
	b1 = quantile(samples, [0.25, 0.75])
	b2 = hpdi(samples, alpha=0.5)

	p1 = plot_density_interval(samples, b1;
	  xlab="Proportion water (p)", title="50% PI");
	p2 = plot_density_interval(samples, b2;
	  xlab="Proportion water (p)", title="50% HPDI");
end;

# ╔═╡ 15bf932c-f787-11ea-295b-6f68b4c100ec
plot(p1, p2, layout=(1, 2))

# ╔═╡ 15cd1a88-f787-11ea-33fa-4f5d7cebab01
md"## End of Fig3.3t.jl"

# ╔═╡ Cell order:
# ╟─bd4d25cc-f786-11ea-18fd-3326d3896012
# ╠═15b6d7c8-f787-11ea-0fd9-4738ba47533d
# ╠═15b712ce-f787-11ea-245f-55620c5d1b05
# ╠═15b79bfe-f787-11ea-010f-5d6a113c0f19
# ╠═15bf2d08-f787-11ea-1316-e71570fded77
# ╠═15bf932c-f787-11ea-295b-6f68b4c100ec
# ╟─15cd1a88-f787-11ea-33fa-4f5d7cebab01
