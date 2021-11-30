### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 5d783968-f503-11ea-2a7f-c178016c9aee
using Pkg, DrWatson

# ╔═╡ 5da06014-f503-11ea-3247-9d56efd51315
using StatisticalRethinking

# ╔═╡ 1d31f6f8-f503-11ea-3be1-471a142ee605
md"## Clip-03-17-19t.jl"

# ╔═╡ 5d9ff58c-f503-11ea-1e45-85013f83f4b2
@quickactivate "StatisticalRethinkingTuring"

# ╔═╡ 5dad7252-f503-11ea-0fcd-2d550e087c4b
md"### snippet 3.11"

# ╔═╡ 5dae18ec-f503-11ea-3d5a-07f5924669ae
begin
	p_grid = range(0, step=0.00001, stop=1)
	prior = ones(length(p_grid))
	likelihood = [pdf(Binomial(3, p), 3) for p in p_grid]
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
	samples = sample(p_grid, Weights(posterior), length(p_grid));
end

# ╔═╡ 5db96b02-f503-11ea-38fe-41d84e68ba6e
md"### snippet 3.14-16"

# ╔═╡ 5dba1d40-f503-11ea-22f2-4d79ce1517fe
begin
	p1 = density(samples;
	  xlab="Proportion water (p)", ylab="Density", lab="density",
	  title="Mean, median & mode", leg=:topleft);
	vline!(p1, [mode(samples)], lab="mode")
	vline!(p1, [median(samples)], lab="median")
	vline!(p1, [mean(samples)], lab="mean")
end

# ╔═╡ 5dc56b70-f503-11ea-0c94-53ba51a38824
md"### snippet 3.17-19"

# ╔═╡ b2bfb106-f503-11ea-2d25-85732a809965
loss = map(p -> sum(posterior .* abs.(p .- p_grid)), p_grid);

# ╔═╡ 5dc62824-f503-11ea-0605-3ddfb6928156
begin
	m = findmin(loss)
	p2 = plot(loss;
	  xlab="Decision (p)", ylab="Expected proportional loss",
	  title="Loss value", lab="Loss function");
	scatter!([m[2]], [m[1]], lab="p_grid[$(m[2])]=$(round(m[1], digits=3))")
end;

# ╔═╡ 5dd1f8e8-f503-11ea-0bfd-759d3c65ed80
plot(p1, p2, layout=(1, 2))

# ╔═╡ 5ddca9c8-f503-11ea-2527-9d51a7cbf3a8
md"## End of clip-03-17-19t.jl"

# ╔═╡ Cell order:
# ╟─1d31f6f8-f503-11ea-3be1-471a142ee605
# ╠═5d783968-f503-11ea-2a7f-c178016c9aee
# ╠═5d9ff58c-f503-11ea-1e45-85013f83f4b2
# ╠═5da06014-f503-11ea-3247-9d56efd51315
# ╟─5dad7252-f503-11ea-0fcd-2d550e087c4b
# ╠═5dae18ec-f503-11ea-3d5a-07f5924669ae
# ╠═5db96b02-f503-11ea-38fe-41d84e68ba6e
# ╠═5dba1d40-f503-11ea-22f2-4d79ce1517fe
# ╟─5dc56b70-f503-11ea-0c94-53ba51a38824
# ╠═b2bfb106-f503-11ea-2d25-85732a809965
# ╠═5dc62824-f503-11ea-0605-3ddfb6928156
# ╠═5dd1f8e8-f503-11ea-0bfd-759d3c65ed80
# ╟─5ddca9c8-f503-11ea-2527-9d51a7cbf3a8
