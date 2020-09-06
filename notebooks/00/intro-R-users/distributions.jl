### A Pluto.jl notebook ###
# v0.11.12

using Markdown
using InteractiveUtils

# ╔═╡ 998407c2-ed46-11ea-1c81-e33bc9e340f7
using Pkg, DrWatson

# ╔═╡ c3112cc0-ed46-11ea-1ecd-b3d8b8d51ac1
using StatisticalRethinking

# ╔═╡ 115d3284-ed4c-11ea-366e-33869314ff98
md"## Distributions.jl"

# ╔═╡ c30a6ffc-ed46-11ea-1036-57f7da453ef2
@quickactivate "StatisticalRethinkingTuring"

# ╔═╡ c3119548-ed46-11ea-34e2-efd188f0457c
# snippet 2.1
begin
	ways = [0, 3, 8, 9, 0]
	ways = ways / sum(ways)
end

# ╔═╡ edb843be-ed46-11ea-1f85-6b7a7bff28ec
md"##### Working with distributions is a little different between Julia and R.

Instead of having `rbinom`, `dbinom` and `pbinom` we just
have the `Binomial` distribution which to work with:"

# ╔═╡ c31a57f2-ed46-11ea-1303-c3bf7b3b7c8f
d = Binomial(9, 0.5)

# ╔═╡ a0b518c8-ed51-11ea-06b1-cd86fe8d1229
md"##### Parameters of Binomial."

# ╔═╡ 5f0e19e6-ed51-11ea-3130-ebd24a5f62bc
fieldnames(Binomial)

# ╔═╡ 3be503f6-ed51-11ea-0958-5fa25106d3f5
md"##### Number of trials parameter."

# ╔═╡ c3272656-ed46-11ea-3971-a3cc353218cb
d.n

# ╔═╡ c872e37e-ed47-11ea-141c-f7a43395df4d
md"##### Singe random draw."

# ╔═╡ c22522ca-ed47-11ea-1fb8-a7e582a8969e
rand(d)

# ╔═╡ 5ff49ae0-ed4c-11ea-0994-ebf5f6c17922
md"##### 9 draws."

# ╔═╡ c3298d88-ed46-11ea-1a28-578706bf1724
rand(d, 9)                          # 

# ╔═╡ c335d4ee-ed46-11ea-3c3b-637d1a3d5bc8
pdf(d, 6)                        # probability density of getting a 6

# ╔═╡ c3403d4e-ed46-11ea-3dce-7906235096c6
cdf(d, 6)                        # cumulative probability of getting 6

# ╔═╡ 17eee314-ed87-11ea-0343-fb0148a66a04
pkg"st"

# ╔═╡ c354cad4-ed46-11ea-1244-9d97d407103f
md"## End of distributions.jl"

# ╔═╡ Cell order:
# ╟─115d3284-ed4c-11ea-366e-33869314ff98
# ╠═998407c2-ed46-11ea-1c81-e33bc9e340f7
# ╠═c30a6ffc-ed46-11ea-1036-57f7da453ef2
# ╠═c3112cc0-ed46-11ea-1ecd-b3d8b8d51ac1
# ╠═c3119548-ed46-11ea-34e2-efd188f0457c
# ╟─edb843be-ed46-11ea-1f85-6b7a7bff28ec
# ╠═c31a57f2-ed46-11ea-1303-c3bf7b3b7c8f
# ╟─a0b518c8-ed51-11ea-06b1-cd86fe8d1229
# ╠═5f0e19e6-ed51-11ea-3130-ebd24a5f62bc
# ╟─3be503f6-ed51-11ea-0958-5fa25106d3f5
# ╠═c3272656-ed46-11ea-3971-a3cc353218cb
# ╟─c872e37e-ed47-11ea-141c-f7a43395df4d
# ╠═c22522ca-ed47-11ea-1fb8-a7e582a8969e
# ╟─5ff49ae0-ed4c-11ea-0994-ebf5f6c17922
# ╠═c3298d88-ed46-11ea-1a28-578706bf1724
# ╠═c335d4ee-ed46-11ea-3c3b-637d1a3d5bc8
# ╠═c3403d4e-ed46-11ea-3dce-7906235096c6
# ╠═17eee314-ed87-11ea-0343-fb0148a66a04
# ╟─c354cad4-ed46-11ea-1244-9d97d407103f
