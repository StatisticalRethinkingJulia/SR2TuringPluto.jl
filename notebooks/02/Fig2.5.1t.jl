### A Pluto.jl notebook ###
# v0.15.1

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
	using Turing
	using StatisticalRethinking
	using PlutoUI
end

# ╔═╡ cb3b80f4-df47-11ea-18e1-dd4b1f2b5cde
md"## Fig 2.5.1t"

# ╔═╡ f65a77d8-df47-11ea-271a-41999fd773fb
md"""
##### This clip is only intended to generate part of Fig 2.5 using a PlutoUI slider.

It is not intended to show how to use Turing (yet)!

This notebook demonstrates simple PlutoUI interactivity."""

# ╔═╡ 147b737a-df48-11ea-3679-77200acb11f0
md"### 1. Create a Turing model:"

# ╔═╡ 0b3fbb40-df48-11ea-08f2-479bc2292d46
@model m2_0t(n, k) = begin
  theta ~ Beta(1, 1) # prior
  k ~ Binomial(n, theta) # model
  return k, theta
end;

# ╔═╡ c07942f0-ec64-11ea-0002-e734a075766d
md"### 2. Generate observed data."

# ╔═╡ 2f43c3b0-df48-11ea-2d13-99adeddbe90a
md"##### n can go from 1:9"

# ╔═╡ 150bbb7c-e0d0-11ea-3ba0-57135fa7c974
@bind n Slider(1:18, default=9)

# ╔═╡ fe925dc4-ec64-11ea-3d14-192e171af40b
md"### 3. Sample for varying n:"

# ╔═╡ 48d028ea-dfcb-11ea-018b-25399925cdef
begin
	dfa2_0t = DataFrame()
	k = sum([1,0,1,1,1,0,1,0,1,1,0,1,1,1,0,1,0,1][1:n])
	chns = sample(m2_0t(n, k), NUTS(0.65), 1000)
	dfa2_0t.theta = reshape(chns[:theta].data, length(chns))
	Text(precis(dfa2_0t; io=String))
end

# ╔═╡ 67e8eb60-ec65-11ea-3e6e-63e8df353e85
md"### 6. Show the posterior."

# ╔═╡ 991b1400-df48-11ea-02a7-e9bdf39427ff
begin
  plot(xlims=(0.0, 1.0), ylims=(0.0, 4.0), leg=false)
  hline!([1.0], line=(:dash))
  density!(dfa2_0t.theta, line=(:dash))
 end

# ╔═╡ a730740e-04ee-11eb-29e5-0381bea25aea
md"##### A few ways to summarize these samples."

# ╔═╡ b95103b2-04ee-11eb-0006-35fdefe0a6d9
md"###### Normal distribution."

# ╔═╡ 8005d0a0-04db-11eb-3d05-35b7c41fee0f
part2_0t = Particles(dfa2_0t)

# ╔═╡ c9213c9e-04ee-11eb-3587-c99b444ec5fa
md"###### MAP based on samples."

# ╔═╡ 4929df5e-04e0-11eb-2b59-f9551aacedea
quap2_0 = quap(dfa2_0t)

# ╔═╡ e1717f20-04ee-11eb-3e82-21ad9de15caa
md"###### Turing quap estimate."

# ╔═╡ e35ad042-04db-11eb-10bc-7f9a46553139
quap2_0t = quap(m2_0t(n, k))

# ╔═╡ 13851f4a-dfc8-11ea-0933-cb4f026bcf42
md"## End of Fig2.5.1t.jl"

# ╔═╡ Cell order:
# ╟─cb3b80f4-df47-11ea-18e1-dd4b1f2b5cde
# ╠═cd2366c0-e0c5-11ea-287a-abc0804397c8
# ╠═9fb491f0-df47-11ea-3cf9-6fa3cee85c33
# ╟─f65a77d8-df47-11ea-271a-41999fd773fb
# ╟─147b737a-df48-11ea-3679-77200acb11f0
# ╠═0b3fbb40-df48-11ea-08f2-479bc2292d46
# ╟─c07942f0-ec64-11ea-0002-e734a075766d
# ╟─2f43c3b0-df48-11ea-2d13-99adeddbe90a
# ╠═150bbb7c-e0d0-11ea-3ba0-57135fa7c974
# ╟─fe925dc4-ec64-11ea-3d14-192e171af40b
# ╠═48d028ea-dfcb-11ea-018b-25399925cdef
# ╟─67e8eb60-ec65-11ea-3e6e-63e8df353e85
# ╠═991b1400-df48-11ea-02a7-e9bdf39427ff
# ╟─a730740e-04ee-11eb-29e5-0381bea25aea
# ╟─b95103b2-04ee-11eb-0006-35fdefe0a6d9
# ╠═8005d0a0-04db-11eb-3d05-35b7c41fee0f
# ╟─c9213c9e-04ee-11eb-3587-c99b444ec5fa
# ╠═4929df5e-04e0-11eb-2b59-f9551aacedea
# ╟─e1717f20-04ee-11eb-3e82-21ad9de15caa
# ╠═e35ad042-04db-11eb-10bc-7f9a46553139
# ╟─13851f4a-dfc8-11ea-0933-cb4f026bcf42
