### A Pluto.jl notebook ###
# v0.11.11

using Markdown
using InteractiveUtils

# ╔═╡ 4133e4a2-ed4d-11ea-11cd-0da48c09c2e8
using Pkg, DrWatson

# ╔═╡ 54bd3938-ed4d-11ea-29e4-c1c7020c1400
begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

# ╔═╡ 0c5fe59e-ed4f-11ea-35a8-45c3aff83362
md"## Broadcasting.jl"

# ╔═╡ 54ccfc74-ed4d-11ea-366d-b738d54ff79a
p_grid = range(0, 1, length = 5)

# ╔═╡ 54d9bce8-ed4d-11ea-279e-3bc3bf465e28
prior = ones(5)

# ╔═╡ 54da91d4-ed4d-11ea-04f5-e9a4f8c659a3
likelihood = pdf.(Binomial.(9, p_grid), 6)

# ╔═╡ 54e45ff4-ed4d-11ea-0288-e3475c709007
posterior = likelihood .* prior

# ╔═╡ 54e8f906-ed4d-11ea-2fee-316c2df1d8f6
posterior ./= sum(posterior)

# ╔═╡ 54f1592c-ed4d-11ea-0210-0be7741d2604
md"##### Broadcasting is another important difference between Julia and R.

You don't have automatic
broadcasting (making functions work the same whether you run it on a single
number or on an array of numbers). However, fixing this is very easy,
you only have to annotate your function call with a dot:
  `f(single_number)` -> `f.(array_of_numbers)`.

That means the translation of the line `likelihood = pdf.(Binomial.(9, p_grid), 6)` goes something like:
1. First or every value in `p_grid` make a binomial distribution with that p value.
2. For every distribution then take the pdf of 6.

The same in the next line, this is elementwise multiplication.

For more infos read [this](https://julialang.org/blog/2017/01/moredots)."

# ╔═╡ e822d97e-ed4f-11ea-2c18-17ecfc81fdbd
md"## End of broadcasting.jl"

# ╔═╡ Cell order:
# ╟─0c5fe59e-ed4f-11ea-35a8-45c3aff83362
# ╠═4133e4a2-ed4d-11ea-11cd-0da48c09c2e8
# ╠═54bd3938-ed4d-11ea-29e4-c1c7020c1400
# ╠═54ccfc74-ed4d-11ea-366d-b738d54ff79a
# ╠═54d9bce8-ed4d-11ea-279e-3bc3bf465e28
# ╠═54da91d4-ed4d-11ea-04f5-e9a4f8c659a3
# ╠═54e45ff4-ed4d-11ea-0288-e3475c709007
# ╠═54e8f906-ed4d-11ea-2fee-316c2df1d8f6
# ╠═54f1592c-ed4d-11ea-0210-0be7741d2604
# ╟─e822d97e-ed4f-11ea-2c18-17ecfc81fdbd
