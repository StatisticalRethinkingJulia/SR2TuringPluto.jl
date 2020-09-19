### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ f252364e-f2e3-11ea-09b5-ff1a685abd56
using Pkg, DrWatson

# ╔═╡ c126b36a-f2e3-11ea-1b0f-a7d38ecffcd8
md"## Clip-03-01t.jl"

# ╔═╡ f2526f60-f2e3-11ea-2d32-276d0b91d917
@quickactivate "StatisticalRethinkingTuring"

# ╔═╡ f252ea76-f2e3-11ea-093a-01c585ee2ea1
md"### snippet 3.1"

# ╔═╡ f25a9eec-f2e3-11ea-1611-9129817fe7a6
begin
	Pr_Positive_Vampire = 0.95
	Pr_Positive_Mortal = 0.01
	Pr_Vampire = 0.001
	Pr_Positive = Pr_Positive_Vampire * Pr_Vampire + Pr_Positive_Mortal * (1 - Pr_Vampire)
	Pr_Vampire_Positive = Pr_Positive_Vampire * Pr_Vampire / Pr_Positive
	Pr_Vampire_Positive
end

# ╔═╡ f25b24a2-f2e3-11ea-302d-f3e7f49f645b
md"## End of clip-03-01t.jl"

# ╔═╡ Cell order:
# ╟─c126b36a-f2e3-11ea-1b0f-a7d38ecffcd8
# ╠═f252364e-f2e3-11ea-09b5-ff1a685abd56
# ╠═f2526f60-f2e3-11ea-2d32-276d0b91d917
# ╟─f252ea76-f2e3-11ea-093a-01c585ee2ea1
# ╠═f25a9eec-f2e3-11ea-1611-9129817fe7a6
# ╟─f25b24a2-f2e3-11ea-302d-f3e7f49f645b
