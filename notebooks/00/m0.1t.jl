### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 8a4bdc28-f9b7-11ea-068e-eb70e542a539
using Pkg, DrWatson

# ╔═╡ 8a4c20a2-f9b7-11ea-2d67-a1e2531dcbd6
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ 4eb6e25c-f9b7-11ea-3e43-05d4d6f124f5
md"## m0.1t.jl"

# ╔═╡ 8a4c985e-f9b7-11ea-22cb-717b854951d6
md"##### Define a simple Normal model with unknown mean and variance."

# ╔═╡ 8a53fb6a-f9b7-11ea-1a81-5b84e8ec7e71
@model gdemo(x, y) = begin
  s ~ InverseGamma(2, 3)
  m ~ Normal(0, sqrt(s))
  x ~ Normal(m, sqrt(s))
  y ~ Normal(m, sqrt(s))
end

# ╔═╡ 8a5473b8-f9b7-11ea-3b00-75cd16c17551
md"#####  Run sampler, collect results."

# ╔═╡ 8a5e1b36-f9b7-11ea-36ef-ebe92bc31e40
chns = mapreduce(c -> sample(gdemo(1.5, 2), NUTS(0.65), 2000), chainscat, 1:4)

# ╔═╡ 8a5ea05e-f9b7-11ea-2cbe-e15b867548b3
md"##### Plot the chains."

# ╔═╡ 8a6d9e26-f9b7-11ea-1d0f-f95263880beb
plot(chns; seriestype=:traceplot)

# ╔═╡ 45425d72-f9b8-11ea-1faa-41e93626b8b5
plot(chns; seriestype=:density)

# ╔═╡ 8a6e37c6-f9b7-11ea-2abd-5db63b061a7d
md"## End m0.1t.jl"

# ╔═╡ Cell order:
# ╟─4eb6e25c-f9b7-11ea-3e43-05d4d6f124f5
# ╠═8a4bdc28-f9b7-11ea-068e-eb70e542a539
# ╠═8a4c20a2-f9b7-11ea-2d67-a1e2531dcbd6
# ╟─8a4c985e-f9b7-11ea-22cb-717b854951d6
# ╠═8a53fb6a-f9b7-11ea-1a81-5b84e8ec7e71
# ╟─8a5473b8-f9b7-11ea-3b00-75cd16c17551
# ╠═8a5e1b36-f9b7-11ea-36ef-ebe92bc31e40
# ╟─8a5ea05e-f9b7-11ea-2cbe-e15b867548b3
# ╠═8a6d9e26-f9b7-11ea-1d0f-f95263880beb
# ╠═45425d72-f9b8-11ea-1faa-41e93626b8b5
# ╟─8a6e37c6-f9b7-11ea-2abd-5db63b061a7d
