# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.8.1
#   kernelspec:
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# +
using CSV
using DataFrames
using Turing
using Logging
using StatisticalRethinking
using StatisticalRethinkingPlots
using StatsBase
using Random
using LaTeXStrings
using StatsPlots

using Dagitty

default(label=false)
Logging.disable_logging(Logging.Warn);
# -

# # 5.1 Spurious association

# Code 5.1

d = DataFrame(CSV.File("data/WaffleDivorce.csv"))
d[!,:D] = standardize(ZScoreTransform, d.Divorce)
d[!,:M] = standardize(ZScoreTransform, d.Marriage)
d[!,:A] = standardize(ZScoreTransform, d.MedianAgeMarriage);

# Code 5.2

std(d.MedianAgeMarriage)

# Code 5.3

# +
Random.seed!(100)

@model function model_m5_1(A, D)
    σ ~ Exponential(1)
    a ~ Normal(0, 0.2)
    bA ~ Normal(0, 0.5)
    μ = @. a + bA * A
    D ~ MvNormal(μ, σ)
end

m5_1 = sample(model_m5_1(d.A, d.D), NUTS(), 1000)
m5_1_df = DataFrame(m5_1)
prior = sample(model_m5_1([0], [0]), Prior(), 1000)
prior_df = DataFrame(prior);
# -

# Code 5.4

# +
# calculate μ for every prior sample on age=-2 and age=2
bounds = [-2, 2]
μ = StatisticalRethinking.link(prior_df, [:a, :bA], bounds)
μ = hcat(μ...);

p = plot(xlab="Median age marriage (std)", ylab="Divorce rate (std)")
for μₚ ∈ first(eachrow(μ), 50)
    plot!(bounds, μₚ; c=:black, alpha=0.3)
end
display(p)
# -

# Code 5.5

# +
A_seq = range(-3, 3.2; length=30)

μ = StatisticalRethinking.link(m5_1_df, [:a, :bA], A_seq)
μ = hcat(μ...)
μ_mean = mean.(eachcol(μ))
μ_PI = PI.(eachcol(μ))
μ_PI = vcat(μ_PI'...)

@df d scatter(:A, :D; xlab="Median age marriage (std)", ylab="Divorce rate (std)")
plot!(A_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.2)
# -

# Code 5.6

# +
Random.seed!(100)

@model function model_m5_2(M, D)
    σ ~ Exponential(1)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    μ = @. a + bM * M
    D ~ MvNormal(μ, σ)
end

m5_2 = sample(model_m5_2(d.M, d.D), NUTS(), 1000)
m5_2_df = DataFrame(m5_2);

# +
M_seq = range(-1.74, 2.8; length=30)

μ = StatisticalRethinking.link(m5_2_df, [:a, :bM], M_seq)
μ = hcat(μ...)
μ_mean = mean.(eachcol(μ))
μ_PI = PI.(eachcol(μ))
μ_PI = vcat(μ_PI'...)

@df d scatter(:M, :D; xlab="Marriage rate (std)", ylab="Divorce rate (std)")
plot!(M_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.2)
# -

# Code 5.7

g = Dagitty.DAG(:A => :M, :A => :D, :M => :D)
drawdag(g, [0, 1, 2], [0, 1, 0])

# Code 5.8

g = Dagitty.DAG(:A => :M, :A => :D)
implied_conditional_independencies(g)

# Code 5.9

g = Dagitty.DAG(:A => :M, :A => :D, :M => :D)
implied_conditional_independencies(g)

# Code 5.10

# +
@model function model_m5_3(A, M, D)
    σ ~ Exponential(1)
    a ~ Normal(0, 0.2)
    bA ~ Normal(0, 0.5)
    bM ~ Normal(0, 0.5)
    μ = @. a + bA * A + bM * M
    D ~ MvNormal(μ, σ)
end

m5_3 = sample(model_m5_3(d.A, d.M, d.D), NUTS(), 1000)
m5_3_df = DataFrame(m5_3)
precis(m5_3_df)
# -

# Code 5.11

coeftab_plot(m5_1_df, m5_2_df, m5_3_df; pars=(:bA, :bM), names=["m5.1", "m5.2", "m5.3"])

# Code 5.12

Random.seed!(100)
N = 50
age = rand(Normal(), N)
mar = rand.(Normal.(-age))
div = rand.(Normal.(age));

s1 = DataFrame(sample(model_m5_1(age, div), NUTS(), 1000))
s2 = DataFrame(sample(model_m5_2(mar, div), NUTS(), 1000))
s3 = DataFrame(sample(model_m5_3(age, mar, div), NUTS(), 1000));
coeftab_plot(s1, s2, s3; pars=(:bA, :bM), names=["s1", "s2", "s3"])

Random.seed!(100)
N = 50
age = rand(Normal(), N)
mar = rand.(Normal.(-age))
div = rand.(Normal.(age .+ mar));

s1 = DataFrame(sample(model_m5_1(age, div), NUTS(), 1000))
s2 = DataFrame(sample(model_m5_2(mar, div), NUTS(), 1000))
s3 = DataFrame(sample(model_m5_3(age, mar, div), NUTS(), 1000));
coeftab_plot(s1, s2, s3; pars=(:bA, :bM), names=["s1", "s2", "s3"])

# Code 5.13

# +
Random.seed!(100)

@model function model_m5_4(A, M)
    σ ~ Exponential(1)
    a ~ Normal(0, 0.2)
    bAM ~ Normal(0, 0.5)
    μ = @. a + bAM * A
    M ~ MvNormal(μ, σ)
end

m5_4 = sample(model_m5_4(d.A, d.M), NUTS(), 1000)
m5_4_df = DataFrame(m5_4);
# -

# Code 5.14

mu = StatisticalRethinking.link(m5_4_df, [:a, :bAM], d.A);
mu = hcat(mu...)
mu_mean = mean.(eachcol(mu))
mu_resid = mu_mean .- d.M;

# +
# Side-note: how to plot the residuals
# getting yerr - list of 2-tuples with distance to the regression line
yerr = collect(zip(-clamp.(mu_resid, -Inf, -0.0), clamp.(mu_resid, 0, Inf)));

plot(d.A, mu_mean; xlab="Age at marriage (std)", ylab="Marriage rate (std)")
scatter!(d.A, d.M)
scatter!(d.A, d.M; yerr=yerr, markersize=0)
# -

# Code 5.15

# +
# explicit link form before I improved it
mu = [
    @. r.a + r.bA * d.A + r.bM * d.M
    for r ∈ eachrow(m5_3_df)
]

mu = vcat(mu'...)
mu_mean = mean.(eachcol(mu))
mu_PI = PI.(eachcol(mu))
mu_PI = vcat(mu_PI'...);

D_sim = [
    rand(MvNormal((@. r.a + r.bA * d.A + r.bM * d.M), r.σ))
    for r ∈ eachrow(m5_3_df)
]
D_sim = vcat(D_sim'...);
D_PI = PI.(eachcol(D_sim))
D_PI = vcat(D_PI'...);
# -

# Code 5.16

yerr = mu_PI[:,2] .- mu_mean
scatter(d.D, mu_mean; xlab="Observed divorce", ylab="Predicted divorce", yerr=yerr)
plot!(x->x; style=:dash)

# Code 5.17

loc_flags = d.Loc .∈ (["ID", "UT", "RI", "ME"],);
loc_idxes = findall(loc_flags);
anns = [
    (d.D[idx] - 0.1, mu_mean[idx], (d.Loc[idx], 8))
    for idx in loc_idxes
]
annotate!(anns)

# Code 5.18

Random.seed!(100)
N = 100
x_real = rand(Normal(), N)
x_spur = rand.(Normal.(x_real))
y = rand.(Normal.(x_real))
df = DataFrame(:y => y, :x_real => x_real, :x_spur => x_spur);

# Code 5.19

# +
d1 = DataFrame(CSV.File("data/WaffleDivorce.csv"))
d = DataFrame(
    :D => standardize(ZScoreTransform, d1.Divorce),
    :M => standardize(ZScoreTransform, d1.Marriage),
    :A => standardize(ZScoreTransform, d1.MedianAgeMarriage),
);

@model function model_m5_3A(A, M, D)
    # A → D ← M
    σ ~ Exponential(1)
    a ~ Normal(0, 0.2)
    bA ~ Normal(0, 0.5)
    bM ~ Normal(0, 0.5)
    μ = @. a + bA * A + bM * M
    D ~ MvNormal(μ, σ)
    # A → M
    σ_M ~ Exponential(1)
    aM ~ Normal(0, 0.2)
    bAM ~ Normal(0, 0.5)
    μ_M = @. aM + bAM * A
    M ~ MvNormal(μ_M, σ_M)
end

m5_3A = sample(model_m5_3A(d.A, d.M, d.D), NUTS(), 1000)
m5_3A_df = DataFrame(m5_3A)
precis(m5_3A_df)
# -

# Code 5.20

A_seq = range(-2, 2; length=30);

# Code 5.21

# +
s_M, s_D = [], []

for r ∈ eachrow(m5_3A_df)
    M = rand(MvNormal((@. r.aM + r.bAM * A_seq), r.σ_M))
    D = rand(MvNormal((@. r.a + r.bA * A_seq + r.bM * M), r.σ))
    push!(s_M, M)
    push!(s_D, D)
end

s_M = vcat(s_M'...)
s_D = vcat(s_D'...);
# -

# Code 5.22

# +
μ_D = mean.(eachcol(s_D))
PI_D = vcat(PI.(eachcol(s_D))'...)

plot(
    A_seq, [μ_D, μ_D]; 
    fillrange=PI_D, fillalpha=0.2, color=:black,
    xlab="manupulated A", ylab="counterfactual D",
    title="Total counterfactual effect of A on D"
)

# +
μ_M = mean.(eachcol(s_M))
PI_M = vcat(PI.(eachcol(s_M))'...)

plot(
    A_seq, [μ_M, μ_M]; 
    fillrange=PI_M, fillalpha=0.2, color=:black,
    xlab="manupulated A", ylab="counterfactual M",
    title="Total counterfactual effect of A on M"
)
# -

# Code 5.23

# +
sim2_A = @. ([20, 30] - 26.1) / 1.24;
s2_M, s2_D = [], []

for r ∈ eachrow(m5_3A_df)
    M = rand(MvNormal((@. r.aM + r.bAM * sim2_A), r.σ_M))
    D = rand(MvNormal((@. r.a + r.bA * sim2_A + r.bM * M), r.σ))
    push!(s2_M, M)
    push!(s2_D, D)
end

s2_M = vcat(s2_M'...)
s2_D = vcat(s2_D'...);
mean(s2_D[:,2] - s2_D[:,1])
# -

# Code 5.24

# +
M_seq = range(-2, 2; length=30)
s_D = []

for r ∈ eachrow(m5_3A_df)
    # A is zero, so, we drop it from the μ term
    D = rand(MvNormal((@. r.a + r.bM * M_seq), r.σ))
    push!(s_D, D)
end

s_D = vcat(s_D'...);

μ_D = mean.(eachcol(s_D))
PI_D = vcat(PI.(eachcol(s_D))'...)

plot(
    M_seq, [μ_D, μ_D]; 
    fillrange=PI_D, fillalpha=0.2, color=:black,
    xlab="manupulated M", ylab="counterfactual D",
    title="Total counterfactual effect of M on D"
)
# -

# Code 5.25

A_seq = range(-2, 2; length=30);

# Code 5.26

s_M = [
    rand(MvNormal((@. r.aM + r.bAM * A_seq), r.σ_M))
    for r ∈ eachrow(m5_3A_df)
]
s_M = vcat(s_M'...);

# Code 5.27

s_D = [
    rand(MvNormal((@. r.a + r.bA * A_seq + r.bM * M), r.σ))
    for (r, M) ∈ zip(eachrow(m5_3A_df), eachrow(s_M))
]
s_D = vcat(s_D'...);

# # 5.2 Masked relationship

# Code 5.28

# +
d = DataFrame(CSV.File("data/milk.csv",  missingstring="NA"))

# get rid of dots in column names
rename!(n -> replace(n, "." => "_"), d)

describe(d)
# -

# Code 5.29

# +
d[!,:K] = standardize(ZScoreTransform, d.kcal_per_g)
d[!,:M] = standardize(ZScoreTransform, log.(d.mass))

# column contains missing values, need to propagate them on standartization
d[!,:N] = d.neocortex_perc
non_miss = findall(!ismissing, d.N);
d[non_miss,:N] = standardize(ZScoreTransform, disallowmissing(d.N[non_miss]));
# -

# Code 5.30

# +
@model function model_m5_5_draft(N, K)
    a ~ Normal(0, 1)
    bN ~ Normal(0, 1)
    σ ~ Exponential(1)
    μ = @. a + bN * N
    K ~ MvNormal(μ, σ)
end

try
    m5_5_draft = sample(model_m5_5_draft(d.N, d.K), NUTS(), 1000)
catch e
    if isa(e, MethodError)
        s = sprint(showerror, e)
        println(s)
    end
end
# -

# Code 5.31

d.neocortex_perc

# Code 5.32

dcc = d[completecases(d[!,[:K,:N,:M]]),:];

# Code 5.33

m5_5_draft = sample(model_m5_5_draft(dcc.N, dcc.K), NUTS(), 1000);

# Code 5.34

# +
prior = sample(model_m5_5_draft(dcc.N, dcc.K), Prior(), 1000)
prior_df = DataFrame(prior)
xseq = [-2, 2]
μ = StatisticalRethinking.link(prior_df, [:a, :bN], xseq)
μ = hcat(μ...);

p = plot(; xlim=xseq, ylim=xseq, 
    xlab="neocortex percent (std)", ylab="kilocal per g (std)", 
    title=L"a \sim \mathcal{N}(0,1), bN \sim \mathcal{N}(0,1)"
)
for y ∈ first(eachrow(μ), 50)
    plot!(p, xseq, y; c=:black, alpha=0.3)
end
p
# -

# Code 5.35

# +
@model function model_m5_5(N, K)
    a ~ Normal(0, 0.2)
    bN ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = @. a + bN * N
    K ~ MvNormal(μ, σ)
end

m5_5 = sample(model_m5_5(dcc.N, dcc.K), NUTS(), 1000)
m5_5_df = DataFrame(m5_5);

# +
prior = sample(model_m5_5(dcc.N, dcc.K), Prior(), 1000)
prior_df = DataFrame(prior)

μ = StatisticalRethinking.link(prior_df, [:a, :bN], xseq)
μ = hcat(μ...);

p2 = plot(; xlim=xseq, ylim=xseq, 
    xlab="neocortex percent (std)", ylab="kilocal per g (std)", 
    title=L"a \sim \mathcal{N}(0,0.2), bN \sim \mathcal{N}(0,0.5)"
)
for y ∈ first(eachrow(μ), 50)
    plot!(p2, xseq, y; c=:black, alpha=0.3)
end
p2
# -

# Code 5.36

precis(m5_5_df)

# Code 5.37

# +
xseq = range(minimum(dcc.N) - 0.15, maximum(dcc.N) + 0.15; length=30)
μ = StatisticalRethinking.link(m5_5_df, [:a, :bN], xseq);
μ = hcat(μ...)
μ_mean = mean.(eachcol(μ))
μ_PI = PI.(eachcol(μ))
μ_PI = vcat(μ_PI'...)

@df dcc scatter(:N, :K; xlab="neocortex percent (std)", ylab="kilocal per g (std)")
plot!(xseq, [μ_mean, μ_mean]; lw=2, fillrange=μ_PI, fillalpha=0.2, color=:black)
# -

# Code 5.38

# +
@model function model_m5_6(M, K)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = @. a + bM * M
    K ~ MvNormal(μ, σ)
end

m5_6 = sample(model_m5_6(dcc.M, dcc.K), NUTS(), 1000)
m5_6_df = DataFrame(m5_6)
precis(m5_6_df)

# +
xseq = range(minimum(dcc.M) - 0.15, maximum(dcc.M) + 0.15; length=30)
μ = StatisticalRethinking.link(m5_6_df, [:a, :bM], xseq);
μ = hcat(μ...)
μ_mean = mean.(eachcol(μ))
μ_PI = PI.(eachcol(μ))
μ_PI = vcat(μ_PI'...)

@df dcc scatter(:M, :K; xlab="log body mass (std)", ylab="kilocal per g (std)")
plot!(xseq, [μ_mean, μ_mean]; lw=2, fillrange=μ_PI, fillalpha=0.2, color=:black)

# +
@model function model_m5_7(N, M, K)
    a ~ Normal(0, 0.2)
    bN ~ Normal(0, 0.5)
    bM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = @. a + bN * N + bM * M
    K ~ MvNormal(μ, σ)
end

m5_7 = sample(model_m5_7(dcc.N, dcc.M, dcc.K), NUTS(), 1000)
m5_7_df = DataFrame(m5_7)
precis(m5_7_df)
# -

# Code 5.40

coeftab_plot(m5_7_df, m5_6_df, m5_5_df; pars=(:bM, :bN), names=("m5.7", "m5.6", "m5.5"))

# Code 5.41
#
# The code in the book corresponds to the bottom-right figure, which keeps N=0 (despite stated in the text).
#
# Below is the code to produce the bottom-left figure (M=0).

# +
xseq = range(minimum(dcc.N) - 0.15, maximum(dcc.N) + 0.15; length=30)
μ = StatisticalRethinking.link(m5_7_df, [:a, :bN], xseq);
μ = hcat(μ...)
μ_mean = mean.(eachcol(μ))
μ_PI = PI.(eachcol(μ))
μ_PI = vcat(μ_PI'...)

plot(title="Counterfactual holding M=0", 
    xlab="neocortex percent (std)", ylab="kilocal per g (std)")
plot!(xseq, [μ_mean, μ_mean]; lw=2, fillrange=μ_PI, fillalpha=0.2, color=:black)

# +
xseq = range(minimum(dcc.M) - 0.15, maximum(dcc.M) + 0.15; length=30)
μ = StatisticalRethinking.link(m5_7_df, [:a, :bM], xseq);
μ = hcat(μ...)
μ_mean = mean.(eachcol(μ))
μ_PI = PI.(eachcol(μ))
μ_PI = vcat(μ_PI'...)

plot(title="Counterfactual holding N=0", 
    xlab="log body mass (std)", ylab="kilocal per g (std)")
plot!(xseq, [μ_mean, μ_mean]; lw=2, fillrange=μ_PI, fillalpha=0.2, color=:black)
# -

# Code 5.42

# M → K ← N
# M → N
Random.seed!(100)
n = 100
M = rand(Normal(), n)
N = [rand(Normal(μ)) for μ ∈ M]
K = [rand(Normal(μ)) for μ ∈ N .- M] 
d_sim = DataFrame(:K => K, :N => N, :M => M);

s5 = sample(model_m5_5(d_sim.N, d_sim.K), NUTS(), 1000)
s6 = sample(model_m5_6(d_sim.M, d_sim.K), NUTS(), 1000)
s7 = sample(model_m5_7(d_sim.N, d_sim.M, d_sim.K), NUTS(), 1000)
s5_df = DataFrame(s5)
s6_df = DataFrame(s6)
s7_df = DataFrame(s7)
coeftab_plot(s7_df, s6_df, s5_df; pars=(:bM, :bN), names=("s7", "s6", "s5"))

# Code 5.43

# +
Random.seed!(100)

# M → K ← N
# N → M
n = 100
N = rand(Normal(), n)
M = [rand(Normal(μ)) for μ ∈ N]
K = [rand(Normal(μ)) for μ ∈ N .- M] 
d_sim2 = DataFrame(:K => K, :N => N, :M => M);

# M → K ← N
# M ← U → N
n = 100
U = rand(Normal(), n)
N = [rand(Normal(μ)) for μ ∈ U]
M = [rand(Normal(μ)) for μ ∈ U]
K = [rand(Normal(μ)) for μ ∈ N .- M] 
d_sim3 = DataFrame(:K => K, :N => N, :M => M);
# -

s5 = sample(model_m5_5(d_sim2.N, d_sim2.K), NUTS(), 1000)
s6 = sample(model_m5_6(d_sim2.M, d_sim2.K), NUTS(), 1000)
s7 = sample(model_m5_7(d_sim2.N, d_sim2.M, d_sim2.K), NUTS(), 1000)
s5_df = DataFrame(s5)
s6_df = DataFrame(s6)
s7_df = DataFrame(s7)
coeftab_plot(s7_df, s6_df, s5_df; pars=(:bM, :bN), names=("s7", "s6", "s5"))

s5 = sample(model_m5_5(d_sim3.N, d_sim3.K), NUTS(), 1000)
s6 = sample(model_m5_6(d_sim3.M, d_sim3.K), NUTS(), 1000)
s7 = sample(model_m5_7(d_sim3.N, d_sim3.M, d_sim3.K), NUTS(), 1000)
s5_df = DataFrame(s5)
s6_df = DataFrame(s6)
s7_df = DataFrame(s7)
coeftab_plot(s7_df, s6_df, s5_df; pars=(:bM, :bN), names=("s7", "s6", "s5"))

# Code 5.44

# +
dag5_7 = Dagitty.DAG(:M => :K, :N => :K, :M => :N)
drawdag(dag5_7, [1, 0, 1], [0, 0, 1])

# equivalentDAGs is TODO in Dagitty.jl
# -

# # 5.3 Categorical variables

# Code 5.45

d = DataFrame(CSV.File("data/Howell1.csv"))
describe(d)

# Code 5.46

cnt = 10_000
μ_female = rand(Normal(178, 20), cnt)
μ_male = rand(Normal(178, 20), cnt) + rand(Normal(0, 10), cnt)
precis(DataFrame(:μ_female => μ_female, :μ_male => μ_male))

# Code 5.47

d[!,:sex] = ifelse.(d.male .== 1, 2, 1)
describe(d.sex)

# Code 5.48

# +
@model function model_m5_8(sex, height)
    σ ~ Uniform(0, 50)
    a ~ MvNormal([178, 178], 20)
    height ~ MvNormal(a[sex], σ)
end

m5_8 = sample(model_m5_8(d.sex, d.height), NUTS(), 1000)
m5_8_df = DataFrame(m5_8)
precis(m5_8_df)
# -

# Code 5.49

m5_8_df[!,:diff_fm] = m5_8_df[:,"a[1]"] - m5_8_df[:,"a[2]"]
precis(m5_8_df)

# Code 5.50

# +
d = DataFrame(CSV.File("data/milk.csv"))

# get rid of dots in column names
rename!(n -> replace(n, "." => "_"), d)

levels(d.clade)
# -

# Code 5.51

d[!,:clade_id] = indexin(d.clade, levels(d.clade));

# Code 5.52

# +
d[!,:K] = standardize(ZScoreTransform, d.kcal_per_g);
clade_counts = maximum(levels(d.clade_id))

@model function model_m5_9(clade_id, K)
    clade_μ = zeros(clade_counts)
    a ~ MvNormal(clade_μ, 0.5)
    σ ~ Exponential(1)
    K ~ MvNormal(a[clade_id], σ)
end

m5_9 = sample(model_m5_9(d.clade_id, d.K), NUTS(), 1000)
m5_9_df = DataFrame(m5_9)
# get rid of square brackets in a params
rename!(n -> replace(n, r"\[|\]" => ""), m5_9_df)

pars = [:a1, :a2, :a3, :a4]
p_names = map(v -> "$(v[1]): $(v[2])", zip(pars, levels(d.clade)))

coeftab_plot(m5_9_df; pars=pars, pars_names=p_names, xlab="expected kcal (std)")
# -

# Code 5.53

# it took me a while to find a seed which make Slytherin to stand out. 
# So it is just a seed, not the model property
Random.seed!(31)
d[!,:house] = sample(1:4, nrow(d));

# Code 5.54

# +
house_counts = maximum(levels(d.house))

@model function model_m5_10(clade_id, house, K)
    clade_μ = zeros(clade_counts)
    house_μ = zeros(house_counts)
    a ~ MvNormal(clade_μ, 0.5)
    h ~ MvNormal(house_μ, 0.5)
    σ ~ Exponential(1)
    μ = a[clade_id] .+ h[house]
    K ~ MvNormal(μ, σ)
end

m5_10 = sample(model_m5_10(d.clade_id, d.house, d.K), NUTS(), 1000)
m5_10_df = DataFrame(m5_10)
precis(m5_10_df)
# -


