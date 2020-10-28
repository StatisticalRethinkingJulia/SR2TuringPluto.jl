
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Turing intro clip 3_quap.jl"

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
    μ = a .+ b .* weights
    for i in eachindex(heights)
        heights[i] ~ Normal(μ[i], σ)
    end
end

begin
	m4_3t = m4_3(df.weight, df.height)
    sampler=NUTS(0.65)
	nsamples=2000
	nchains=4
    chns4_3t = mapreduce(c -> sample(m4_3t, sampler, nsamples), chainscat, 1:nchains)
end

q4_3t = quap(m4_3t)

begin
	scatter(df.weight, df.height, leg=false)
	quap4_3t = DataFrame(rand(q4_3t.distr, 10_000)', q4_3t.params)
	a_map = mean(quap4_3t.a)
	b_map = mean(quap4_3t.b)
	plot!(x, a_map .+ b_map .* x)
end

