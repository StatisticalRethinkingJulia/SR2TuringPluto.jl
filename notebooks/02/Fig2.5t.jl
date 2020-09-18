### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ cd2366c0-e0c5-11ea-287a-abc0804397c8
using Pkg, DrWatson

# ╔═╡ 9fb491f0-df47-11ea-3cf9-6fa3cee85c33
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ cb3b80f4-df47-11ea-18e1-dd4b1f2b5cde
md"## Fig 2.5t"

# ╔═╡ f65a77d8-df47-11ea-271a-41999fd773fb
md"""
##### This clip is only intended to generate part of Fig 2.5."""

# ╔═╡ 147b737a-df48-11ea-3679-77200acb11f0
md"##### Create a Turing modelt:"

# ╔═╡ 0b3fbb40-df48-11ea-08f2-479bc2292d46
@model globe_toss(W, L) = begin
    p ~ Uniform(0, 1)
    W ~ Binomial(W + L, p)
end

# ╔═╡ c080ab9a-eede-11ea-1a02-7747110ce510
m = globe_toss(6, 3);

# ╔═╡ c0ac2b1c-eede-11ea-0af7-c90c3dc45666
r = quap(m)

# ╔═╡ c0b70690-eede-11ea-1086-95193b22f266
d = Normal(r.coef.p, √collect(reshape(r.vcov, 1, 1))[1])

# ╔═╡ c0c5c68a-eede-11ea-3754-471209a14604
s = rand(d, 10000);

# ╔═╡ c0d26d52-eede-11ea-394a-153cc7000dae
h1 = histogram(s, normalize = :pdf, lab="sample density")

# ╔═╡ 2f325012-f91f-11ea-16d0-034a3f1f9f4b
begin
	k = 3
	n = 6
	chn = sample(globe_toss(n, k), NUTS(0.65), 11000)
end

# ╔═╡ d64cb1f6-f921-11ea-3d8c-cd654a1bd6f6
plot(chn; seriestype=:traceplot)

# ╔═╡ f6bbd200-f921-11ea-3c1d-2b7617b99656
begin
	plot(chn; seriestype=:density)
	histogram!(s, normalize = :pdf, lab="sample density")
end

# ╔═╡ 13851f4a-dfc8-11ea-0933-cb4f026bcf42
md"## End of Fig2.5t.jl"

# ╔═╡ Cell order:
# ╠═cb3b80f4-df47-11ea-18e1-dd4b1f2b5cde
# ╠═cd2366c0-e0c5-11ea-287a-abc0804397c8
# ╠═9fb491f0-df47-11ea-3cf9-6fa3cee85c33
# ╟─f65a77d8-df47-11ea-271a-41999fd773fb
# ╟─147b737a-df48-11ea-3679-77200acb11f0
# ╠═0b3fbb40-df48-11ea-08f2-479bc2292d46
# ╠═c080ab9a-eede-11ea-1a02-7747110ce510
# ╠═c0ac2b1c-eede-11ea-0af7-c90c3dc45666
# ╠═c0b70690-eede-11ea-1086-95193b22f266
# ╠═c0c5c68a-eede-11ea-3754-471209a14604
# ╠═c0d26d52-eede-11ea-394a-153cc7000dae
# ╠═2f325012-f91f-11ea-16d0-034a3f1f9f4b
# ╠═d64cb1f6-f921-11ea-3d8c-cd654a1bd6f6
# ╠═f6bbd200-f921-11ea-3c1d-2b7617b99656
# ╟─13851f4a-dfc8-11ea-0933-cb4f026bcf42
