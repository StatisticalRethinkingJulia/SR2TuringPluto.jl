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
md"## Pluto-helpers-prettytable.jl"

# ╔═╡ b4656386-2db1-11eb-089f-993d79325858
begin
	# Define the experiment
	niter = 4000
	nparams = 3
	nchains = 2

	# some sample experiment results
	val = randn(niter, nparams, nchains) .+ [1, 2, 3]'
	val = hcat(val, rand(1:2, niter, 1, nchains))

	# construct a Chains object
	chns = Chains(val, start = 1, thin = 2)
end;

# ╔═╡ d1670910-2c0d-11eb-12d0-616c518e22da
Text(sprint(show, "text/plain", chns))

# ╔═╡ 5e493ac4-2db2-11eb-0567-7d635483d547
begin
	bounds = hpd(chns)
	Text(sprint(show, "text/plain", bounds))
end

# ╔═╡ 36733c28-2da0-11eb-3b0b-d5c9b85e16ba
CHNS(chns)

# ╔═╡ eff124b4-2dda-11eb-2717-dde3f4d416e6
HPD(chns)

# ╔═╡ 3ce617ae-f934-11ea-0e4a-ff1ec9ee2610
md"## End pluto-helpers-prettytable.jl"

# ╔═╡ Cell order:
# ╟─83ad95c0-f933-11ea-1552-23ac7681803c
# ╠═3c57ceb0-f934-11ea-1bef-c16b272d0c6b
# ╠═3c580a7e-f934-11ea-3b0c-4d5d7d5a7382
# ╠═b4656386-2db1-11eb-089f-993d79325858
# ╠═d1670910-2c0d-11eb-12d0-616c518e22da
# ╠═5e493ac4-2db2-11eb-0567-7d635483d547
# ╠═36733c28-2da0-11eb-3b0b-d5c9b85e16ba
# ╠═eff124b4-2dda-11eb-2717-dde3f4d416e6
# ╟─3ce617ae-f934-11ea-0e4a-ff1ec9ee2610
