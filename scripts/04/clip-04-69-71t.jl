
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## clip-04-69-71t"


md"### snippet 4.64"

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	scale!(df, [:weight])
end;

md"### snippet 4.69"

f_cube(weight_s, a, b1, b2, b3) = a + b1 * weight_s + b2 * weight_s^2 + b3 * weight_s^3

@model function cube(weight_s, heights)
    a ~ Normal(178, 20)
    b1 ~ LogNormal(0, 1)
    b2 ~ Normal(0, 10)
    b3 ~ Normal(0, 10)
    σ ~ Uniform(0, 50)
    μ = f_cube.(weight_s, a, b1, b2, b3)
    heights ~ MvNormal(μ, σ)
end

begin
	quap4_6t = quap(cube(df.weight_s, df.height), NelderMead())
	quap4_6t.coef
end

md"### snippets 4.70, 4.71"

begin
	weight_seq_2 = range(-2.2, 2, length = 30)
	post_2 = DataFrame(rand(quap4_6t.distr, 1_000)', quap4_6t.params)
	mu_2 = f_cube.(weight_seq_2, post_2.a', post_2.b1', post_2.b2', post_2.b3') |> meanlowerupper
	sim_2 = rand.(Normal.(mu_2.raw, post_2.σ')) |> meanlowerupper
	weight_seq_rescaled = weight_seq_2 .* std(df.weight) .+ mean(df.weight)
	scatter(df.weight, df.height, ms = 3, alpha = 0.7, legend = false)
	plot!(weight_seq_rescaled, mu_2.mean, ribbon = (mu_2.mean .- mu_2.lower, mu_2.upper .- mu_2.mean))
	plot!(weight_seq_rescaled, sim_2.lower, fillrange = sim_2.upper, alpha = 0.3, la = 0.0, c = 2)
end

md"## End of clip-04-69-71t.jl"

