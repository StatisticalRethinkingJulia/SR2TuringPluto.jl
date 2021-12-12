### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ eb033006-f788-11ea-086f-357e2353381f
using Pkg, DrWatson

# ╔═╡ eb03721e-f788-11ea-10b4-07a4a4032d76
begin
	@quickactivate "SR2TuringPluto"
	using StatisticalRethinking
end

# ╔═╡ c4b79432-f788-11ea-1092-df50f9dca4bf
md"## Fig3.5t.jl"

# ╔═╡ eb03e064-f788-11ea-274c-cbe1c1978108
begin
  N = 100000
  d = rand(Binomial(9, 0.7), N);
  histogram(d; normalize=:probability, 
    bins=-0.5:1:9.5, leg=false, xticks=0:9, bar_width=0.2)
end

# ╔═╡ eb0aecd8-f788-11ea-14da-fdfc08b9a986
md"## End of Fig3.5t.jl"

# ╔═╡ Cell order:
# ╟─c4b79432-f788-11ea-1092-df50f9dca4bf
# ╠═eb033006-f788-11ea-086f-357e2353381f
# ╠═eb03721e-f788-11ea-10b4-07a4a4032d76
# ╠═eb03e064-f788-11ea-274c-cbe1c1978108
# ╟─eb0aecd8-f788-11ea-14da-fdfc08b9a986
