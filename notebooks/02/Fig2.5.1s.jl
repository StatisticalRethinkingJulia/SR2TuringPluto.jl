### A Pluto.jl notebook ###
# v0.11.12

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ cd2366c0-e0c5-11ea-287a-abc0804397c8
using Pkg, DrWatson

# ╔═╡ 9fb491f0-df47-11ea-3cf9-6fa3cee85c33
begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
	using Turing
	using PlutoUI
	include(srcdir() * "/quap.jl")
end

# ╔═╡ cb3b80f4-df47-11ea-18e1-dd4b1f2b5cde
md"## Fig 2.5.1s"

# ╔═╡ f65a77d8-df47-11ea-271a-41999fd773fb
md"""
##### This clip is only intended to generate part of Fig 2.5 using a PlutoUI slider.

It is not intended to show how to use Stan (yet)!"""

# ╔═╡ 147b737a-df48-11ea-3679-77200acb11f0
md"##### Create a Turing modelt:"

# ╔═╡ 0b3fbb40-df48-11ea-08f2-479bc2292d46
@model globethrowing(W, L) = begin
    p ~ Uniform(0, 1)
    W ~ Binomial(W + L, p)
end

# ╔═╡ c080ab9a-eede-11ea-1a02-7747110ce510
m = globethrowing(6, 3)

# ╔═╡ c0ac2b1c-eede-11ea-0af7-c90c3dc45666
r = quap(m)

# ╔═╡ c0b70690-eede-11ea-1086-95193b22f266
d = MvNormal([r.coef.p], collect(reshape(r.vcov, 1, 1)))

# ╔═╡ c0c5c68a-eede-11ea-3754-471209a14604
s = rand(d, 10000)'

# ╔═╡ c0d26d52-eede-11ea-394a-153cc7000dae
histogram(collect(s), normalize = :probability)

# ╔═╡ c0df0566-eede-11ea-0141-99e33e557afd
r.vcov[1] |> sqrt

# ╔═╡ c0f56444-eede-11ea-216c-152682bc8e97
histogram(rand(d, 10000)', normalize = :probability)

# ╔═╡ 2331e85c-df48-11ea-1551-b54d9e48188c


# ╔═╡ c07942f0-ec64-11ea-0002-e734a075766d
md"### 2. Generate observed data."

# ╔═╡ 2f43c3b0-df48-11ea-2d13-99adeddbe90a
md"##### n can go from 1:9"

# ╔═╡ 150bbb7c-e0d0-11ea-3ba0-57135fa7c974
@bind n Slider(1:18, default=9)

# ╔═╡ fe925dc4-ec64-11ea-3d14-192e171af40b
md"### 3. Create a stan_sample data object (a Dict):"

# ╔═╡ 48d028ea-dfcb-11ea-018b-25399925cdef
begin
	k = [1,0,1,1,1,0,1,0,1,1,0,1,1,1,0,1,0,1][1:n]
  	m2_0_data = Dict("n" => n, "k" => sum(k[1:n]));
end

# ╔═╡ 32d38400-ec65-11ea-22cc-d1d82d7d4b08
md"### 4. Sample posterior."

# ╔═╡ 4d624f04-ec65-11ea-2640-d5a76b39e6ba
md"### 5. If successful, retieve the draws."

# ╔═╡ 67e8eb60-ec65-11ea-3e6e-63e8df353e85
md"### 6. Show the posterior."

# ╔═╡ 991b1400-df48-11ea-02a7-e9bdf39427ff
begin
  plot(xlims=(0.0, 1.0), ylims=(0.0, 4.0), leg=false)
  hline!([1.0], line=(:dash))
  density!(dfs.theta, line=(:dash))
 end

# ╔═╡ 13851f4a-dfc8-11ea-0933-cb4f026bcf42
md"## End of Fig2.5.1s.jl"

# ╔═╡ Cell order:
# ╟─cb3b80f4-df47-11ea-18e1-dd4b1f2b5cde
# ╠═cd2366c0-e0c5-11ea-287a-abc0804397c8
# ╠═9fb491f0-df47-11ea-3cf9-6fa3cee85c33
# ╟─f65a77d8-df47-11ea-271a-41999fd773fb
# ╠═147b737a-df48-11ea-3679-77200acb11f0
# ╠═0b3fbb40-df48-11ea-08f2-479bc2292d46
# ╠═c080ab9a-eede-11ea-1a02-7747110ce510
# ╠═c0ac2b1c-eede-11ea-0af7-c90c3dc45666
# ╠═c0b70690-eede-11ea-1086-95193b22f266
# ╠═c0c5c68a-eede-11ea-3754-471209a14604
# ╠═c0d26d52-eede-11ea-394a-153cc7000dae
# ╠═c0df0566-eede-11ea-0141-99e33e557afd
# ╠═c0f56444-eede-11ea-216c-152682bc8e97
# ╟─2331e85c-df48-11ea-1551-b54d9e48188c
# ╟─c07942f0-ec64-11ea-0002-e734a075766d
# ╟─2f43c3b0-df48-11ea-2d13-99adeddbe90a
# ╠═150bbb7c-e0d0-11ea-3ba0-57135fa7c974
# ╟─fe925dc4-ec64-11ea-3d14-192e171af40b
# ╠═48d028ea-dfcb-11ea-018b-25399925cdef
# ╟─32d38400-ec65-11ea-22cc-d1d82d7d4b08
# ╟─4d624f04-ec65-11ea-2640-d5a76b39e6ba
# ╟─67e8eb60-ec65-11ea-3e6e-63e8df353e85
# ╠═991b1400-df48-11ea-02a7-e9bdf39427ff
# ╟─13851f4a-dfc8-11ea-0933-cb4f026bcf42
