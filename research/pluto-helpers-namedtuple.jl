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
md"## Pluto-helpers-namedtuple.jl"

# ╔═╡ 3c9b1db4-f934-11ea-3d97-f51c2acc9037
md"## snippet 2.6"

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
	Text(precis(quap2_0t_df; io=String))
end

# ╔═╡ d1670910-2c0d-11eb-12d0-616c518e22da
begin
	post2_0t = sample(m2_0t, NUTS(), 5000)
	Text(sprint(show, "text/plain", post2_0t))
end

# ╔═╡ 0604d19a-2cd4-11eb-3797-975be539c612
QM(q2_0t)

# ╔═╡ 57af25e0-2d97-11eb-0c3e-bb88e690e5ae
q2_0t

# ╔═╡ 3ce617ae-f934-11ea-0e4a-ff1ec9ee2610
md"## End pluto_helpers-namedtuple.jl"

# ╔═╡ Cell order:
# ╟─83ad95c0-f933-11ea-1552-23ac7681803c
# ╠═3c57ceb0-f934-11ea-1bef-c16b272d0c6b
# ╠═3c580a7e-f934-11ea-3b0c-4d5d7d5a7382
# ╟─3c9b1db4-f934-11ea-3d97-f51c2acc9037
# ╠═3ca3ac90-f934-11ea-1a8e-edfa2c664f1a
# ╠═3caf45dc-f934-11ea-2d1d-796707d89526
# ╠═3cb18732-f934-11ea-3ebf-991c7b602d12
# ╠═3cbf6b30-f934-11ea-3b10-433589e36648
# ╠═d1670910-2c0d-11eb-12d0-616c518e22da
# ╠═0604d19a-2cd4-11eb-3797-975be539c612
# ╠═57af25e0-2d97-11eb-0c3e-bb88e690e5ae
# ╟─3ce617ae-f934-11ea-0e4a-ff1ec9ee2610
