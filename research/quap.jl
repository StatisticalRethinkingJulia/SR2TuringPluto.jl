using Pkg, DrWatson

@quickactivate "SR2TuringPluto"
using Turing
using StatisticalRethinking

Turing.setprogress!(false)

df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
df = df[df.age .>= 18, :]

@model function ppl4_1(heights)
    μ ~ Normal(178, 20)
    σ ~ Uniform(0, 50)
    heights .~ Normal(μ, σ)
end

m4_1t = ppl4_1(df.height)
precis(m4_1t)

res = quap(m4_1t)
