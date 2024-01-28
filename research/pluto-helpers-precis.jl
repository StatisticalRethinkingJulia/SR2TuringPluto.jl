### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ 3c57ceb0-f934-11ea-1bef-c16b272d0c6b
using Pkg, DrWatson

# ╔═╡ 3c580a7e-f934-11ea-3b0c-4d5d7d5a7382
begin
	@quickactivate "SR2TuringPluto"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ 83ad95c0-f933-11ea-1552-23ac7681803c
md"## pluto-helpers-precis.jl"

# ╔═╡ 3ca3ac90-f934-11ea-1a8e-edfa2c664f1a
@model function ppl2_0(W, L)
    p ~ Uniform(0, 1)
    W ~ Binomial(W + L, p)
end

# ╔═╡ 3caf45dc-f934-11ea-2d1d-796707d89526
m2_0t = ppl2_0(6, 3);

# ╔═╡ 3cb18732-f934-11ea-3ebf-991c7b602d12
q2_0t = quap(m2_0t)

# ╔═╡ 3cbf6b30-f934-11ea-3b10-433589e36648
begin
	quap2_0t_df = DataFrame(:μ => rand(q2_0t.distr, 4000))
	res = Text(precis(quap2_0t_df; io=String))
end

# ╔═╡ 36733c28-2da0-11eb-3b0b-d5c9b85e16ba
describe(quap2_0t_df)

# ╔═╡ 3ce617ae-f934-11ea-0e4a-ff1ec9ee2610
md"## End pluto-helpers-precis.jl"

# ╔═╡ Cell order:
# ╟─83ad95c0-f933-11ea-1552-23ac7681803c
# ╠═3c57ceb0-f934-11ea-1bef-c16b272d0c6b
# ╠═3c580a7e-f934-11ea-3b0c-4d5d7d5a7382
# ╠═3ca3ac90-f934-11ea-1a8e-edfa2c664f1a
# ╠═3caf45dc-f934-11ea-2d1d-796707d89526
# ╠═3cb18732-f934-11ea-3ebf-991c7b602d12
# ╠═3cbf6b30-f934-11ea-3b10-433589e36648
# ╠═36733c28-2da0-11eb-3b0b-d5c9b85e16ba
# ╟─3ce617ae-f934-11ea-0e4a-ff1ec9ee2610
