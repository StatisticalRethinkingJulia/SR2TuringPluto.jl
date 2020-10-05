
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## clip-04-64-79t"


md"### snippet 4.64"

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	scale!(df, [:weight])
end;

f_parabola(weight_s, a, b1, b2) = a + b1 * weight_s + b2 * weight_s^2

@model function parabola(weight_s, heights)
    a ~ Normal(178, 20)
    b1 ~ LogNormal(0, 1)
    b2 ~ Normal(0, 1)
    σ ~ Uniform(0, 50)
    μ = f_parabola.(weight_s, a, b1, b2)
    heights ~ MvNormal(μ, σ)
end

begin
	quap4_5t = quap(parabola(df.weight_s, df.height), NelderMead())
	quap4_5t.coef
end

quap4_5t.vcov

md"### snippet 4.67"

begin
	weight_seq_1 = range(-2.2, 2, length = 30)
	post_1 = DataFrame(rand(quap4_5t.distr, 1_000)', quap4_5t.params)
	mu_1 = f_parabola.(weight_seq_1, post_1.a', post_1.b1', post_1.b2')
	sim_1 = rand.(Normal.(mu_1, post_1.σ'))
	mu_1 = meanlowerupper(mu_1)
	sim_1 = meanlowerupper(sim_1)
end;

md"### snippet 4.68"

begin
	scatter(df.weight_s, df.height, ms = 3, alpha = 0.7, legend = false)
	plot!(weight_seq_1, mu_1.mean, ribbon = (mu_1.mean .- mu_1.lower, mu_1.upper .- mu_1.mean))
	plot!(weight_seq_1, sim_1.lower, fillrange = sim_1.upper, alpha = 0.3, linealpha = 0.0, c = 2)
end

md"## End of clip-04-64-68t.jl"

