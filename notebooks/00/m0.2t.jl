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
	using StatsPlots
end

# ╔═╡ a296146e-f9b8-11ea-255b-0d1b0515811c
using MLJ, MLJModels

# ╔═╡ 4eb6e25c-f9b7-11ea-3e43-05d4d6f124f5
md"## m0.2t.jl"

# ╔═╡ 9de23082-04c7-11eb-3b7c-393a11a68cf4
md"##### This notebook takes a look at the [rstar()](https://arxiv.org/pdf/2003.07900.pdf) diagnostic. Currently I can't get MLJ & MLJModels to work (serialization issue?)."

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

# ╔═╡ 6afed64c-04c8-11eb-1a4f-75b21ce91da9
md"##### Below will likely fail!"

# ╔═╡ d049cfcc-f9b8-11ea-1ed4-3fc6fb56241c
chn ... # sampling results of multiple chains

# ╔═╡ d04a0ae6-f9b8-11ea-3053-21a1d21bb400
md"##### Select classifier used to compute the diagnostic."

# ╔═╡ d04f0cee-f9b8-11ea-354e-e973a232f153
classif = @load XGBoostClassifier

# ╔═╡ d05641e4-f9b8-11ea-01bd-c9c850b9b2c3
md"##### stimate diagnostic."

# ╔═╡ d065bd40-f9b8-11ea-2546-3bc0c2eb0fcd
Rs = rstar(classif, chn)

# ╔═╡ d0682076-f9b8-11ea-30bf-61044f5b34ca
R = mean(Rs)

# ╔═╡ d075a64c-f9b8-11ea-1d7e-0b226c48ba95
md"##### Visualize distribution."

# ╔═╡ d07800ae-f9b8-11ea-045b-53bd1cb7f748
histogram(Rs)

# ╔═╡ 8a6e37c6-f9b7-11ea-2abd-5db63b061a7d
md"## End m0.2t.jl"

# ╔═╡ Cell order:
# ╟─4eb6e25c-f9b7-11ea-3e43-05d4d6f124f5
# ╟─9de23082-04c7-11eb-3b7c-393a11a68cf4
# ╠═8a4bdc28-f9b7-11ea-068e-eb70e542a539
# ╠═8a4c20a2-f9b7-11ea-2d67-a1e2531dcbd6
# ╟─8a4c985e-f9b7-11ea-22cb-717b854951d6
# ╠═8a53fb6a-f9b7-11ea-1a81-5b84e8ec7e71
# ╟─8a5473b8-f9b7-11ea-3b00-75cd16c17551
# ╠═8a5e1b36-f9b7-11ea-36ef-ebe92bc31e40
# ╟─8a5ea05e-f9b7-11ea-2cbe-e15b867548b3
# ╠═8a6d9e26-f9b7-11ea-1d0f-f95263880beb
# ╠═45425d72-f9b8-11ea-1faa-41e93626b8b5
# ╟─6afed64c-04c8-11eb-1a4f-75b21ce91da9
# ╠═a296146e-f9b8-11ea-255b-0d1b0515811c
# ╠═d049cfcc-f9b8-11ea-1ed4-3fc6fb56241c
# ╠═d04a0ae6-f9b8-11ea-3053-21a1d21bb400
# ╠═d04f0cee-f9b8-11ea-354e-e973a232f153
# ╠═d05641e4-f9b8-11ea-01bd-c9c850b9b2c3
# ╠═d065bd40-f9b8-11ea-2546-3bc0c2eb0fcd
# ╠═d0682076-f9b8-11ea-30bf-61044f5b34ca
# ╠═d075a64c-f9b8-11ea-1d7e-0b226c48ba95
# ╠═d07800ae-f9b8-11ea-045b-53bd1cb7f748
# ╟─8a6e37c6-f9b7-11ea-2abd-5db63b061a7d
