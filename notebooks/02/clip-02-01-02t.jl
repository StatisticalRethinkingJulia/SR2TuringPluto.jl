### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ 206337dc-f928-11ea-04ac-278c8fc8079d
using Pkg, DrWatson

# ╔═╡ 20646558-f928-11ea-3351-bf12703651b2
begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

# ╔═╡ 0e0d2086-f926-11ea-0470-451b528c9d02
md"## Clip-02-01-02t.jl"

# ╔═╡ 2069098c-f928-11ea-1008-d75cba2186b4
md"## snippet 2.1"

# ╔═╡ 20699ba4-f928-11ea-3f37-8359c913281d
begin
	ways = [0, 3, 8, 9, 0]
	ways = ways / sum(ways)
end

# ╔═╡ 20762c66-f928-11ea-32ec-53742fbaff9b
md"## snippet 2.2"

# ╔═╡ 207be6ea-f928-11ea-32df-cd2ed1f344d1
md"##### With Distributions.jl working with distributions is a little different from R. Instead of having `rbinom`, `dbinom` and `pbinom` we just got the `Binomial` distribution which to work with."

# ╔═╡ aa0f1104-f928-11ea-0830-cbd5ac5bcb63
d = Binomial(9, 0.5)             # Binomial distribution

# ╔═╡ aa0f4d9a-f928-11ea-038e-c31458a61022
d.n                              # Number of trials parameter

# ╔═╡ aa0fd1cc-f928-11ea-1da8-c5f3559d2d0a
rand(d)                          # singe random draw

# ╔═╡ aa1b33a8-f928-11ea-328e-bf25c307896b
rand(d, 10)                      # 10 random draws

# ╔═╡ aa1bc3c0-f928-11ea-2dff-950e37ca54d3
pdf(d, 6)                        # probability density of getting a 6

# ╔═╡ aa261070-f928-11ea-3502-83a01308f7bd
cdf(d, 6)                        # cumulative probability of getting 6

# ╔═╡ aa271b96-f928-11ea-162b-390eddab0c11
pdf(Binomial(9, 0.5), 6)         # probability density of getting a 6

# ╔═╡ 9e7664be-f928-11ea-20b2-653cf2037050
md"## End of clip-02-01-02t.jl"

# ╔═╡ Cell order:
# ╟─0e0d2086-f926-11ea-0470-451b528c9d02
# ╠═206337dc-f928-11ea-04ac-278c8fc8079d
# ╠═20646558-f928-11ea-3351-bf12703651b2
# ╟─2069098c-f928-11ea-1008-d75cba2186b4
# ╠═20699ba4-f928-11ea-3f37-8359c913281d
# ╟─20762c66-f928-11ea-32ec-53742fbaff9b
# ╟─207be6ea-f928-11ea-32df-cd2ed1f344d1
# ╠═aa0f1104-f928-11ea-0830-cbd5ac5bcb63
# ╠═aa0f4d9a-f928-11ea-038e-c31458a61022
# ╠═aa0fd1cc-f928-11ea-1da8-c5f3559d2d0a
# ╠═aa1b33a8-f928-11ea-328e-bf25c307896b
# ╠═aa1bc3c0-f928-11ea-2dff-950e37ca54d3
# ╠═aa261070-f928-11ea-3502-83a01308f7bd
# ╠═aa271b96-f928-11ea-162b-390eddab0c11
# ╟─9e7664be-f928-11ea-20b2-653cf2037050
