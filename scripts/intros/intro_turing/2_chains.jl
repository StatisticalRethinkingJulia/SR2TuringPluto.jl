
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing 
	using StatisticalRethinking
	Turing.turnprogress(false);
end;

md"## Intro clip 2_chains.jl"

begin
	N = 200
	df = DataFrame()
	df.weight = rand(Uniform(20, 90), N)
	f(x) = mean(df.weight) + 1.6x + rand(Normal(0, 10))
	df.height = f.(df.weight)
	Text(precis(df; io=String))
end

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

plot(chns4_3t; seriestype=:traceplot)

plot(chns4_3t; seriestype=:density)

md"## End of Intro clip 2_chains.jl"

