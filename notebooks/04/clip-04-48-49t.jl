### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 2b78a920-070d-11eb-194e-bd44cc52fd06
using Pkg, DrWatson

# ╔═╡ 0fc3846a-070e-11eb-1144-4918bcd157f3
begin
	@quickactivate "SR2TuringPluto"
	using Turing
	using StatisticalRethinking
end;

# ╔═╡ dfe57196-070c-11eb-1a81-fdf2e39f38ae
md"## Clip-04-48-49t.jl"

# ╔═╡ 2d8739d4-070d-11eb-1bd2-c9a5fd864e4d
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)     # either
end;

# ╔═╡ 37566494-070d-11eb-0c3b-694670d5f2ea
@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* (weights .- x̄)
    heights ~ MvNormal(μ, σ)
end

# ╔═╡ d987cc24-070c-11eb-12ba-4b04b7f9c218
md"### snippets 4.48, 4.49"

# ╔═╡ 04032cd8-070d-11eb-3610-0b5e5c16c8ec
begin
	figs = Vector{Plots.Plot{Plots.GRBackend}}(undef, 3)
		for (indx, N) in enumerate([50, 150, 352])

		dN = df[1:N, :]
		quapN = quap(m4_3(dN.weight, dN.height))

		post = rand(quapN.distr, 20)
		post = DataFrame(post', quapN.params)

		figs[indx] = scatter(dN.weight, dN.height)
		for p in eachrow(post)
			plot!(x, p.a .+ p.b .* (x .- mean(dN.weight)), color = "black", alpha = 0.3)
		end
		plot!(legend = false, xlabel = "weight", ylabel = "height")
	end
	plot(figs..., layout=(1, 3))
end

# ╔═╡ b8f032bc-070f-11eb-2ad9-c33312592d49
md"## End of clip-04-48-49t.jl"

# ╔═╡ Cell order:
# ╟─dfe57196-070c-11eb-1a81-fdf2e39f38ae
# ╠═2b78a920-070d-11eb-194e-bd44cc52fd06
# ╠═0fc3846a-070e-11eb-1144-4918bcd157f3
# ╠═2d8739d4-070d-11eb-1bd2-c9a5fd864e4d
# ╠═37566494-070d-11eb-0c3b-694670d5f2ea
# ╟─d987cc24-070c-11eb-12ba-4b04b7f9c218
# ╠═04032cd8-070d-11eb-3610-0b5e5c16c8ec
# ╟─b8f032bc-070f-11eb-2ad9-c33312592d49
