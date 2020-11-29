# m4_3t.jl

using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking
Turing.setprogress!(false)

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)
end;

@model function ppl4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* (weights .- x̄)
    for i in eachindex(heights)
        heights[i] ~ Normal(μ[i], σ)
    end
end

m4_3t = ppl4_3(df.weight, df.height)
q4_3t = quap(m4_3t, NelderMead())

nchains = 4; sampler = NUTS(0.65); nsamples=2000
chns4_3t = mapreduce(c -> sample(m4_3t, sampler, nsamples), chainscat, 1:nchains)
chns4_3t_df = DataFrame(chns4_3t)

quap4_3t_df = DataFrame(rand(q4_3t.distr, 10_000)', q4_3t.params)

part4_3t =	Particles(chns4_3t[[:a, :b, :σ]])

# End of m4.3t.jl
