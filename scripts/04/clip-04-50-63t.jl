
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Clip-04-50-63t.jl"

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)     # either
end;

@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* (weights .- x̄)
    heights ~ MvNormal(μ, σ)
end

m4_3t = m4_3(df.weight, df.height)
quap4_3t = quap(m4_3t)

md"### snippets 4.50 - 4.52"

begin
	post = DataFrame(rand(quap4_3t.distr, 1_000)', quap4_3t.params)
	mu_at_50 = post.a + post.b * (50 - x̄)
end;

density(mu_at_50)

quantile(mu_at_50, (0.1, 0.9))

md"### snippet 4.53 - 4.55"

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

md"### snippets 4.59 - 4.61"

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

md"### snippet 4.62"

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

md"## End of clip-04-50-63t.jl"

