### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 38253cd6-f36c-11ea-397f-3d3dce9bac78
using Pkg, DrWatson

# ╔═╡ 38257d6a-f36c-11ea-3d0a-c19b131da11c
begin
	@quickactivate "SR2TuringPluto"
	using StatisticalRethinking
end

# ╔═╡ 80e5d0da-f36b-11ea-19fa-d9c7e930c6b1
md"## Clip-03-11-16t.jl"

# ╔═╡ 3825f68c-f36c-11ea-07a7-470cc146cc6e
md"### snippet 3.11"

# ╔═╡ 382ee7ec-f36c-11ea-198a-eb0a137d43d1
begin
	p_grid = range(0, step=0.001, stop=1)
	prior = ones(length(p_grid))
	likelihood = [pdf(Binomial(3, p), 3) for p in p_grid]
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
end;

# ╔═╡ 382f603a-f36c-11ea-1f05-7716e241d580
md"##### Draw 10000 samples from this posterior distribution."

# ╔═╡ 383b80ce-f36c-11ea-23a0-97b799c3315a
begin
	N = 10000
	samples = sample(p_grid, Weights(posterior), N);
end;

# ╔═╡ 383c1926-f36c-11ea-3f27-5562a255ce30
md"### snippet 3.13"

# ╔═╡ 38455cac-f36c-11ea-3634-61af92f66f5c
hpdi(samples, alpha=0.11)

# ╔═╡ 3846bb8a-f36c-11ea-086f-23705940b3d4
md"### snippet 3.14"

# ╔═╡ 38517d70-f36c-11ea-1dd8-83f037833053
mode(samples)

# ╔═╡ 3852341a-f36c-11ea-20ed-75fbf044d1d2
md"### snippet 3.15"

# ╔═╡ 38581a4c-f36c-11ea-381e-8dc2daabea5a
mean(samples)

# ╔═╡ 3863c502-f36c-11ea-1f72-ed1eb6c862fb
md"### snippet 3.16"

# ╔═╡ 386500de-f36c-11ea-307d-f3d81413fa54
median(samples)

# ╔═╡ 386c6f04-f36c-11ea-38e1-ebc6cdda48ef
md"##### Plot density."

# ╔═╡ 387396e4-f36c-11ea-3f58-ed848b00749a
begin
	density(samples, lab="density")
	vline!(hpdi(samples, alpha=0.5), line=:dash, lab="hpdi")
	vline!(quantile(samples, [0.25, 0.75]), line=:dash, lab="quantile (pi)")
end

# ╔═╡ 387a79fa-f36c-11ea-3c91-379b9c7b4ab3
md"## End of clip-03-11-16t.jl"

# ╔═╡ Cell order:
# ╟─80e5d0da-f36b-11ea-19fa-d9c7e930c6b1
# ╠═38253cd6-f36c-11ea-397f-3d3dce9bac78
# ╠═38257d6a-f36c-11ea-3d0a-c19b131da11c
# ╟─3825f68c-f36c-11ea-07a7-470cc146cc6e
# ╠═382ee7ec-f36c-11ea-198a-eb0a137d43d1
# ╟─382f603a-f36c-11ea-1f05-7716e241d580
# ╠═383b80ce-f36c-11ea-23a0-97b799c3315a
# ╟─383c1926-f36c-11ea-3f27-5562a255ce30
# ╠═38455cac-f36c-11ea-3634-61af92f66f5c
# ╟─3846bb8a-f36c-11ea-086f-23705940b3d4
# ╠═38517d70-f36c-11ea-1dd8-83f037833053
# ╟─3852341a-f36c-11ea-20ed-75fbf044d1d2
# ╠═38581a4c-f36c-11ea-381e-8dc2daabea5a
# ╟─3863c502-f36c-11ea-1f72-ed1eb6c862fb
# ╠═386500de-f36c-11ea-307d-f3d81413fa54
# ╟─386c6f04-f36c-11ea-38e1-ebc6cdda48ef
# ╠═387396e4-f36c-11ea-3f58-ed848b00749a
# ╟─387a79fa-f36c-11ea-3c91-379b9c7b4ab3
