using DrWatson
@quickactivate "StatReth"

# %%
using DataFrames
using CSV
# using StatsBase
using Distributions
using Turing

include(srcdir("quap.jl"))

# %% 4.26
d = DataFrame(CSV.File(datadir("exp_raw/Howell_1.csv")))
d2 = d[d.age .>= 18, :]

# %% 4.27 - 4.29
@model function height(heights)
    μ ~ Normal(178, 20)
    σ ~ Uniform(0, 50)
    heights .~ Normal(μ, σ)
end

m4_1 = height(d2.height)
q4_1 = quap(m4_1)

# precis(m4_1)

# %% 4.30
start = [mean(Normal(178, 20)), mean(Uniform(0, 50))]
# currently broken https://github.com/TuringLang/Turing.jl/issues/1298
# m4_1 = quap(m, start)

# %% 4.31
@model function height(heights)
    μ ~ Normal(178, 0.1)
    σ ~ Uniform(0, 50)
    heights .~ Normal(μ, σ)
end

m4_2 = height(d2.height)
q4_2 = quap(m4_2, NelderMead())
# precis(m4_2)

# %% 4.32, 4.33
q4_1.vcov

diag(q4_1.vcov)
cov2cor(Matrix(q4_1.vcov), sqrt.(diag(q4_1.vcov)))

# %% 4.34 - 4.36
post = rand(q4_1.distr, 10_000)
post = DataFrame(post', ["μ", "σ"])                         # either
post = DataFrame(rand(q4_1.distr, 10_000)', q4_1.params)    # or

# precis(post)

post = DataFrame(rand(q4_1.distr, 10_000)', ["μ", "σ"])
