using DrWatson
@quickactivate "StatReth"

# %%
using DataFrames
using CSV
using Turing
using Plots

include(srcdir("quap.jl"))
include(srcdir("tools.jl"))

# %% 5.28, 5.29
d = DataFrame!(CSV.File(datadir("exp_raw/milk.csv"), missingstring = "NA"))
# We get the "confusing error message" even earlier because zscore doesn't know how to work
# with missing data, so we'll drop them right out the gate.
dropmissing!(d, [:kcal_per_g, :neocortex_perc, :mass])

d.K = zscore(d.kcal_per_g)
d.N = zscore(d.neocortex_perc)
d.M = zscore(log.(d.mass))

# %% 5.30
# This is why this *would* work on the first try, so we are skipping this

# %% 5.31, 5.32
d.neocortex_perc
# Already did this so it's not going to do anything additional here.
dropmissing!(d, [:kcal_per_g, :neocortex_perc, :mass])

# %% 5.33
@model function milk1(N, K)
    a ~ Normal(0, 1)
    bN ~ Normal(0, 1)
    σ ~ Exponential(1)
    μ = lin(a, bN, N)
    K ~ MvNormal(μ, σ)
end

m5_5_draft = milk1(d.N, d.K)
q5_5_draft = quap(m5_5_draft)

# %% 5.34
prior = sample(m5_5_draft, Prior(), 50) |> DataFrame
xseq = -2:0.1:2
μ = lin(prior.a', prior.bN', xseq)

plot(legend = false, ylims = extrema(xseq))
for c in eachcol(μ)
    plot!(xseq, c, color = :black, alpha = 0.3)
end
plot!()

# %% 5.35
@model function milk2(N, K)
    a ~ Normal(0, 0.2)
    bN ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, bN, N)
    K ~ MvNormal(μ, σ)
end

m5_5 = milk2(d.N, d.K)
q5_5 = quap(m5_5)

# %%
prior = sample(m5_5, Prior(), 50) |> DataFrame
xseq = -2:0.1:2
μ = lin(prior.a', prior.bN', xseq)

plot(legend = false, ylims = extrema(xseq))
for c in eachcol(μ)
    plot!(xseq, c, color = :black, alpha = 0.3)
end
plot!()

# %% 5.36
DataFrame(rand(q5_5.distr, 1000)', q5_5.params) |> precis

# %% 5.37
xseq = range(minimum(d.N) - 0.15, maximum(d.N) + 0.15, length = 30)
post = DataFrame(rand(q5_5.distr, 1000)', q5_5.params)
μ = lin(post.a', post.bN', xseq) |> meanlowerupper

scatter(d.N, d.K, legend = false)
plot!(xseq, μ.mean, ribbon = (μ.mean .- μ.lower, μ.upper .- μ.mean))

# %% 5.38
@model function milk3(M, K)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, bM, M)
    K ~ MvNormal(μ, σ)
end

m5_6 = milk3(d.M, d.K)
q5_6 = quap(m5_6)
post = DataFrame(rand(q5_6.distr, 1000)', q5_6.params)
precis(post)

# %% 5.39
@model function milk4(M, N, K)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    bN ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, bM, M, bN, N)
    K ~ MvNormal(μ, σ)
end

q5_7 = quap(milk4(d.M, d.N, d.K))
post = DataFrame(rand(q5_7.distr, 1000)', q5_7.params)
precis(post)

# %% 5.40
# TODO coefplot

# %% 5.41
xseq = range(minimum(d.M) - 0.15, maximum(d.M) + 0.15, length = 30)
mu = lin(post.a', post.bM', xseq, post.bN', zeros(length(xseq))) |> meanlowerupper
plot(xseq, mu.mean, ribbon = (mu.upper .- mu.mean, mu.mean .- mu.lower))

# %% 5.42
# M -> K <- N
# M -> N
n = 100
M = randn(n)
N = rand.(Normal.(M))
K = rand.(Normal.(N .- M))
d_sim = DataFrame((; K, N, M))

# %% 5.43
# M -> K <- N
# N -> M
n = 100
N = randn(n)
M = rand.(Normal.(N))
K = rand.(Normal.(N .- M))
d_sim2 = DataFrame((; K, N, M))

# M -> K <- N
# M -> U <- N
n = 100
U = randn(n)
N = rand.(Normal.(U))
M = rand.(Normal.(U))
K = rand.(Normal.(N .- M))
d_sim2 = DataFrame((; K, N, M))

# %% 5.44
# TODO working with DAGs
