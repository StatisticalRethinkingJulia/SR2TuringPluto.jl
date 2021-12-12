### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 3c84ab1c-f504-11ea-25f0-6f2b3cbf472b
using Pkg, DrWatson

# ╔═╡ 3c84dde2-f504-11ea-1b8c-95f180a567f3
begin
	@quickactivate "SR2TuringPluto"
	using StatisticalRethinking
end

# ╔═╡ f89aef6a-f503-11ea-306e-e983e6719b50
md"## Clip-03-20-24t.jl"

# ╔═╡ 3c85558c-f504-11ea-3eb4-afb3815d4382
md"### snippet 3.20"

# ╔═╡ 3c8ceeee-f504-11ea-1413-eb6c84b6c00e
pdf.(Binomial(2, 0.7), 0:2)

# ╔═╡ 3c906b5a-f504-11ea-3892-835bb1c47540
md"### snippet 3.21"

# ╔═╡ 3ca25e0a-f504-11ea-3c6b-11d90448ccf0
rand(Binomial(2, 0.7))

# ╔═╡ 3ca30f4e-f504-11ea-3c4e-33ee204e71de
md"### snippet 3.22"

# ╔═╡ 3caea962-f504-11ea-2856-3da15f13e8c1
rand(Binomial(2, 0.7), 10)

# ╔═╡ 3caf7ab8-f504-11ea-0d8b-fb1bff09a842
begin
	N = 100000
	d = rand(Binomial(2, 0.7), N);
	[length(filter(e -> e == i, d)) for i in 0:2] / N |> display

	d = rand(Binomial(9, 0.7), N);
	h1 = histogram([filter(e -> e == i, d) for i in 0:9];
	 bins=-0.5:1:9.5, color=:lightblue, leg=false, xticks=0:9, bar_width=0.2)
end

# ╔═╡ 3cbca094-f504-11ea-3dbb-c7d342341fb2
md"## End of clip-03-20-24t.jl"

# ╔═╡ Cell order:
# ╟─f89aef6a-f503-11ea-306e-e983e6719b50
# ╠═3c84ab1c-f504-11ea-25f0-6f2b3cbf472b
# ╠═3c84dde2-f504-11ea-1b8c-95f180a567f3
# ╟─3c85558c-f504-11ea-3eb4-afb3815d4382
# ╠═3c8ceeee-f504-11ea-1413-eb6c84b6c00e
# ╟─3c906b5a-f504-11ea-3892-835bb1c47540
# ╠═3ca25e0a-f504-11ea-3c6b-11d90448ccf0
# ╟─3ca30f4e-f504-11ea-3c4e-33ee204e71de
# ╠═3caea962-f504-11ea-2856-3da15f13e8c1
# ╠═3caf7ab8-f504-11ea-0d8b-fb1bff09a842
# ╟─3cbca094-f504-11ea-3dbb-c7d342341fb2
