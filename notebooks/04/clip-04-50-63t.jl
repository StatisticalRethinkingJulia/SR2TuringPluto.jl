### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ e0573e4e-0710-11eb-0c38-f54ceea8c90b
using Pkg, DrWatson

# ╔═╡ 9060d3be-0710-11eb-1686-196534333f23
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ 8d881ea4-0710-11eb-1eb7-2337f465fcff
md"## Clip-04-50-63t.jl"

# ╔═╡ 8fd7b794-0710-11eb-2b61-fd8f3fa5951f
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)     # either
end;

# ╔═╡ b639f15a-0725-11eb-3b15-c35097e6acd5
@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* (weights .- x̄)
    heights ~ MvNormal(μ, σ)
end

# ╔═╡ dab57768-0778-11eb-3d4b-43d2eff46c65
begin
	m4_3t = m4_3(df.weight, df.height)
	Text(precis(m4_3t; io=String))
end

# ╔═╡ 183d7c0c-0779-11eb-087b-59e6cb72a2a2
begin
	quap4_3t = quap(m4_3t)
	Text(precis(quap4_3t; io=String))
end

# ╔═╡ 8f6a472e-0710-11eb-143d-6dfc45f27694
md"### snippets 4.50 - 4.52"

# ╔═╡ d463b7e4-070f-11eb-3b02-6f0bf21b511b
begin
	post = DataFrame(rand(quap4_3t.distr, 1_000)', quap4_3t.params)
	mu_at_50 = post.a + post.b * (50 - x̄)
end;

# ╔═╡ bb0877a2-0724-11eb-0ab5-3b1e14ee69ec
density(mu_at_50)

# ╔═╡ bb08bb74-0724-11eb-0e0e-fd511a507191
quantile(mu_at_50, (0.1, 0.9))

# ╔═╡ bb093fe6-0724-11eb-2a6a-2fc93e608313
md"### snippet 4.53 - 4.55"

# ╔═╡ 0613f12c-0725-11eb-15c7-abaa84650f38
begin
	weight_seq = 25:70

	# It's a little unfortunate that you have to write out the formula you have already put
	# into the model. I don't have a better way at the moment though.
	mu = post.a' .+ post.b' .* (weight_seq .- x̄)

	scatter(weight_seq, mu[:, 1:100], legend = false, c = 1, alpha = 0.1)

	# %% 4.56, 4.57
	mu_m = mean.(eachrow(mu))
	mu_lower = quantile.(eachrow(mu), 0.055)
	mu_upper = quantile.(eachrow(mu), 0.945)

	scatter(df.weight, df.height, ms = 3)
	plot!(weight_seq, mu_m, ribbon = (mu_m .- mu_lower, mu_upper .- mu_m))

	# or
	#mu = meanlowerupper(mu)
	#scatter(df.weight, df.height, ms = 3)
	#plot!(weight_seq, mu.mean, ribbon = (mu.mean .- mu.lower, mu.upper .- mu.mean))

end

# ╔═╡ 06146c9c-0725-11eb-078b-b7b2c199795a
md"### snippets 4.59 - 4.61"

# ╔═╡ bb1a7ca4-0724-11eb-2259-c34225bb22f3
begin
	mu1 = meanlowerupper(mu)
	# This isn't really pretty either. I'm not sure how to put this into a function since
	# there isn't a way to know how to create `predict_model` from `height`.
	chn = Chains(Array(post), ["a", "b", "σ"])
	predict_model = m4_3(weight_seq, missing)
	sim1 = predict(predict_model, chn)
	sim1 = meanlowerupper(Array(sim1)')

	scatter(df.weight, df.height, ms = 3, legend = false)
	plot!(weight_seq, mu1.mean, ribbon = (mu1.mean .- mu1.lower, mu1.upper .- mu1.mean))
	plot!(weight_seq, sim1.lower, fillrange = sim1.upper, alpha = 0.3, linealpha = 0.0, c = 2)
end

# ╔═╡ bb1b2adc-0724-11eb-3b86-e972178eee8a
md"### snippet 4.62"

# ╔═╡ bb293b7c-0724-11eb-0a07-3f2bb565e951
begin
	sim2 = predict(predict_model, Chains(rand(quap4_3t.distr, 10_000)', ["a", "b", "σ"])) |> Array
	sim2 = meanlowerupper(sim2')

	scatter(df.weight, df.height, ms = 3, legend = false)
	plot!(weight_seq, mu1.mean, ribbon = (mu1.mean .- mu1.lower, mu1.upper .- mu1.mean))
	plot!(weight_seq, sim2.lower, fillrange = sim2.upper, alpha = 0.3, linealpha = 0.0, c = 2)

	# %% 4.63
	post2 = DataFrame(rand(quap4_3t.distr, 1_000)', quap4_3t.params)
	#weight_seq = 25:70
	normals = Normal.(post2.a' .+ post2.b' .* (weight_seq .- x̄), post2.σ')
	sim3 = rand.(normals)
	sim3 = meanlowerupper(sim3)

	scatter(df.weight, df.height, ms = 3, legend = false)
	plot!(weight_seq, mu1.mean, ribbon = (mu1.mean .- mu1.lower, mu1.upper .- mu1.mean))
	plot!(weight_seq, sim3.lower, fillrange = sim3.upper, alpha = 0.3, linealpha = 0.0, c = 2)
end

# ╔═╡ a0eead00-0724-11eb-0bda-2fce139eebef
md"## End of clip-04-50-63t.jl"

# ╔═╡ Cell order:
# ╠═8d881ea4-0710-11eb-1eb7-2337f465fcff
# ╠═e0573e4e-0710-11eb-0c38-f54ceea8c90b
# ╠═9060d3be-0710-11eb-1686-196534333f23
# ╠═8fd7b794-0710-11eb-2b61-fd8f3fa5951f
# ╠═b639f15a-0725-11eb-3b15-c35097e6acd5
# ╠═dab57768-0778-11eb-3d4b-43d2eff46c65
# ╠═183d7c0c-0779-11eb-087b-59e6cb72a2a2
# ╟─8f6a472e-0710-11eb-143d-6dfc45f27694
# ╠═d463b7e4-070f-11eb-3b02-6f0bf21b511b
# ╠═bb0877a2-0724-11eb-0ab5-3b1e14ee69ec
# ╠═bb08bb74-0724-11eb-0e0e-fd511a507191
# ╟─bb093fe6-0724-11eb-2a6a-2fc93e608313
# ╠═0613f12c-0725-11eb-15c7-abaa84650f38
# ╟─06146c9c-0725-11eb-078b-b7b2c199795a
# ╠═bb1a7ca4-0724-11eb-2259-c34225bb22f3
# ╟─bb1b2adc-0724-11eb-3b86-e972178eee8a
# ╠═bb293b7c-0724-11eb-0a07-3f2bb565e951
# ╟─a0eead00-0724-11eb-0bda-2fce139eebef
