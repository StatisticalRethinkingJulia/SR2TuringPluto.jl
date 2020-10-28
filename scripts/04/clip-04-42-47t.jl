
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Clip-04-42-47t.jl"

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)     # either
end;

md"### snippet 4.42"

@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* (weights .- x̄)
    heights ~ MvNormal(μ, σ)
end

m4_3t = m4_3(df.weight, df.height);

q4_3t = quap(m4_3t, NelderMead())

md"### snippet 4.43"

@model function m4_3l(weights, heights)
    a ~ Normal(178, 20)
    log_b ~ Normal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ exp(log_b) .* (weights .- mean(weights))
    heights .~ Normal.(μ, σ)
end

m4_3tl = m4_3l(df.weight, df.height)

q4_3tl = quap(m4_3tl, NelderMead())

md"### snippets 4.44, 4.45"

round.(q4_3t.vcov, digits = 3)

md"### snippet 4.46"

begin
	scatter(df.weight, df.height, leg=false)
	post = DataFrame(rand(q4_3t.distr, 10_000)', q4_3t.params)
	a_map = mean(post.a)
	b_map = mean(post.b)
	plot!(x, a_map .+ b_map .* (x .- x̄))
end

md"### snippets 4.47"

Text(precis(post; io=String))

md"## End of clip-04-42-47t.jl"

