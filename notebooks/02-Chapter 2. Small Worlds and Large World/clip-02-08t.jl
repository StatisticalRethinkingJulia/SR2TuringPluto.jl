### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 7eeb28ae-f2e1-11ea-2933-955cee17a9bd
using Pkg, DrWatson

# ╔═╡ 7eeb666e-f2e1-11ea-1237-db0c2818822c
begin
	using Distributions
	using StatsPlots, Plots
	using Turing
	using Logging
	using LaTeXStrings
end

# ╔═╡ 6688dfb8-f2df-11ea-2511-2d9736b27265
md"## Clip-02-08t.jl"

# ╔═╡ 7eebee5e-f2e1-11ea-237e-89059c0b9245
md"### snippet 2.8"

# ╔═╡ 7ef86274-f2e1-11ea-25be-5f2b1538b966
md"##### Simple Metropolis algorithm."

# ╔═╡ 7ef91c84-f2e1-11ea-1ed4-c3ef76963769
begin
	n_samples = 10000
	a3d = ones(n_samples,1,1)
	w = 6; l = 3; n = w +l
	p = [0.5]
	for i in 2:n_samples
  		p_new = rand(Normal(p[i-1], 0.1), 1)[1]
  		if  p_new < 0
    		p_new = abs(p_new)
  		end
  		if p_new > 1
    		p_new = 2 - p_new
  		end
  		q0 = pdf(Binomial(n, p[i-1]), w)
  		q1 = pdf(Binomial(n, p_new), w)
  		append!(p, [rand(Uniform(0, 1), 1)[1] < q1/q0 ? p_new : p[i-1]])
	end
end

# ╔═╡ 7f05a540-f2e1-11ea-0d4a-f7c33efa0832
md"##### Create an MCMCChains.Chains object. This Chains object has 10000 samples, one variable and a single chain."

# ╔═╡ 7f0651fe-f2e1-11ea-00f7-79322fafaaf9
begin
	a3d[:, 1, 1] = p
	chns = MCMCChains.Chains(a3d, ["p",])
	Text(sprint(show, "text/plain", summarize(chns)))
end

# ╔═╡ 7f0bac38-f2e1-11ea-0817-093aba68676d
md"##### Show density and computed conjugate solution."

# ╔═╡ 7f16b116-f2e1-11ea-14fb-1ffcca62978a
begin
	x = 0:0.01:1
	density(chns, lab="Samples")
	plot!( x, pdf.(Beta( w+1 , n-w+1 ) , x ), lab="Conjugate solution")
end

# ╔═╡ 7f175f44-f2e1-11ea-31bf-0fbdc4545246
md"## End of clip-02-08t.jl"

# ╔═╡ Cell order:
# ╠═6688dfb8-f2df-11ea-2511-2d9736b27265
# ╠═7eeb28ae-f2e1-11ea-2933-955cee17a9bd
# ╠═7eeb666e-f2e1-11ea-1237-db0c2818822c
# ╟─7eebee5e-f2e1-11ea-237e-89059c0b9245
# ╟─7ef86274-f2e1-11ea-25be-5f2b1538b966
# ╠═7ef91c84-f2e1-11ea-1ed4-c3ef76963769
# ╟─7f05a540-f2e1-11ea-0d4a-f7c33efa0832
# ╠═7f0651fe-f2e1-11ea-00f7-79322fafaaf9
# ╟─7f0bac38-f2e1-11ea-0817-093aba68676d
# ╠═7f16b116-f2e1-11ea-14fb-1ffcca62978a
# ╟─7f175f44-f2e1-11ea-31bf-0fbdc4545246
