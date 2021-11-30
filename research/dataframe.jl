using Pkg, DrWatson

@quickactivate "SR2TuringPluto"
using Turing 
using StatisticalRethinking
Turing.setprogress!(false);

using DataFrames, Turing

N = 200
df = DataFrame()
df.weight = rand(Uniform(20, 90), N)
f(x) = mean(df.weight) + 1.6x + rand(Normal(0, 10))
df.height = f.(df.weight)

@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* weights
    for i in eachindex(heights)
        heights[i] ~ Normal(μ[i], σ)
    end
end

m4_3t = m4_3(df.weight, df.height)
chns4_3t = sample(m4_3t, NUTS(100, 0.65), N);
chns4_3t |> display

# Use DataFrame for easy handling

#= From README docs MCMCChains
Similarly, for DataFrames:

DataFrame(chns)                                 <- Iteration, chain index added
DataFrame(chns[:s])                             <- Still true?
DataFrame(chns, [:parameters])                  <- Still true?
...
=#

Array(chns4_3t, [:parameters]) |> display
DataFrame(Chains(chns4_3t, [:parameters]) |> display
