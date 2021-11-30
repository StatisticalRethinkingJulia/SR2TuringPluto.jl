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
using GLM
using CSV
using Random
using StatsBase
using DataFrames
using Dagitty
using Turing
using StatsPlots
using StatisticalRethinking
using StatisticalRethinkingPlots
using Logging

default(labels=false)
Logging.disable_logging(Logging.Warn);
# -

# Code 6.1

# +
Random.seed!(1917)
N = 200   # grant proposals
p = 0.1   # proportion to select

# uncorrelated newsworthiness and trustworthiness
nw = rand(Normal(), N)
tw = rand(Normal(), N)

# select top 10% of combined score
s = nw .+ tw
q = quantile(s, 1-p)
selected = s .>= q
cor(tw[selected], nw[selected])
# -

scatter(nw[.!selected], tw[.!selected]; xlab="newsworthiness", ylab="trustworthiness", label="not selected")
scatter!(nw[selected], tw[selected]; label="selected")

# # 6.1 Multicollinearity

# Code 6.2

Random.seed!(100)
N = 100
height = rand(Normal(10, 2), N)
leg_prop = rand(Uniform(0.4, 0.5), N)
leg_left = leg_prop .* height .+ rand(Normal(0, 0.02), N)
leg_right = leg_prop .* height .+ rand(Normal(0, 0.02), N)
d = DataFrame(:height => height, :leg_left => leg_left, :leg_right => leg_right);

# Code 6.3

# +
@model function model_m6_1(leg_left, leg_right, height)
    a ~ Normal(10, 100)
    bl ~ Normal(2, 10)
    br ~ Normal(2, 10)
    μ = @. a + bl * leg_left + br * leg_right
    σ ~ Exponential(1)
    height ~ MvNormal(μ, σ)
end

m6_1 = sample(model_m6_1(d.leg_left, d.leg_right, d.height), NUTS(), 1000)
m6_1_df = DataFrame(m6_1)
precis(m6_1_df)
# -

# Code 6.4

coeftab_plot(m6_1_df)

# Code 6.5

scatter(m6_1_df.br, m6_1_df.bl; alpha=0.1)

# Code 6.6

@df m6_1_df density(:br + :bl; lw=2, xlab="sum of bl and br")

# Code 6.7

# +
@model function model_m6_2(leg_left, height)
    a ~ Normal(10, 100)
    bl ~ Normal(2, 10)
    μ = @. a + bl * leg_left
    σ ~ Exponential(1)
    height ~ MvNormal(μ, σ)
end

m6_2 = sample(model_m6_2(d.leg_left, d.height), NUTS(), 1000)
m6_2_df = DataFrame(m6_2)
precis(m6_2_df)
# -

std(m6_1_df.bl), std(m6_1_df.br), std(m6_1_df.bl + m6_1_df.br)

# Code 6.8

# +
d = DataFrame(CSV.File("data/milk.csv",  missingstring="NA"))

# get rid of dots in column names
rename!(n -> replace(n, "." => "_"), d)

d[!,:K] = standardize(ZScoreTransform, d.kcal_per_g)
d[!,:F] = standardize(ZScoreTransform, d.perc_fat)
d[!,:L] = standardize(ZScoreTransform, d.perc_lactose);
# -

# Code 6.9

# +
@model function model_m6_3(F, K)
    a ~ Normal(0, 0.2)
    bF ~ Normal(0, 0.5)
    μ = @. a + F * bF
    σ ~ Exponential(1)
    K ~ MvNormal(μ, σ)
end

m6_3 = sample(model_m6_3(d.F, d.K), NUTS(), 1000)
m6_3_df = DataFrame(m6_3)

@model function model_m6_4(L, K)
    a ~ Normal(0, 0.2)
    bL ~ Normal(0, 0.5)
    μ = @. a + L * bL
    σ ~ Exponential(1)
    K ~ MvNormal(μ, σ)
end

m6_4 = sample(model_m6_4(d.L, d.K), NUTS(), 1000)
m6_4_df = DataFrame(m6_4)

precis(m6_3_df)
precis(m6_4_df)
# -

# Code 6.10

# +
@model function model_m6_5(F, L, K)
    a ~ Normal(0, 0.2)
    bF ~ Normal(0, 0.5)
    bL ~ Normal(0, 0.5)
    μ = @. a + F * bF + L * bL
    σ ~ Exponential(1)
    K ~ MvNormal(μ, σ)
end

m6_5 = sample(model_m6_5(d.F, d.L, d.K), NUTS(), 1000)
m6_5_df = DataFrame(m6_5)
precis(m6_5_df)
# -

# Code 6.11

@df d corrplot([:kcal_per_g :perc_fat :perc_lactose]; seriestype=:scatter, bins=10, grid=false)

# Code 6.12

# +
# get mean stderr for linear model's scale
function stderr_for_r(r)
    σ = sqrt(1-r^2)*var(d.perc_fat)
    fat_scaled = r .* d.perc_fat
    stderr_x = [
        begin
            x = d.perc_fat .+ rand(MvNormal(fat_scaled, σ))
            # add the intercept to the model
            X = hcat(ones(length(x)), x)
            m = lm(X, d.kcal_per_g)
            stderror(m)[2]
        end
        for _ in 1:100
    ]
    s = mean(stderr_x)
end

r_seq = range(0, 0.99; step=0.01)
s = stderr_for_r.(r_seq)
plot(r_seq, s; lw=2, xlab="correlation")
# -

# # 6.2 Post-treatment bias

# Code 6.13

# +
Random.seed!(70)
# number of plants
N = 100
h0 = rand(Normal(10, 2), N)
treatment = repeat(0:1, inner=div(N, 2))
fungus = [rand(Binomial(1, 0.5 - treat*0.4)) for treat in treatment]
h1 = h0 .+ rand(MvNormal(5 .- 3 .* fungus, 1))

d = DataFrame(:h0 => h0, :h1 => h1, :treatment => treatment, :fungus => fungus)
precis(d)
# -

# Code 6.14

sim_p = rand(LogNormal(0, 0.25), 10_000)
precis(DataFrame(:sim_p => sim_p))

# Code 6.15

# +
@model function model_m6_6(h0, h1)
    p ~ LogNormal(0, 0.25)
    σ ~ Exponential(1)
    μ = h0 .* p
    h1 ~ MvNormal(μ, σ)
end

m6_6 = sample(model_m6_6(d.h0, d.h1), NUTS(), 1000)
m6_6_df = DataFrame(m6_6)
precis(m6_6_df)
# -

# Code 6.16

# +
@model function model_m6_7(h0, treatment, fungus, h1)
    a ~ LogNormal(0, 0.2)
    bt ~ Normal(0, 0.5)
    bf ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    p = @. a + bt*treatment + bf*fungus
    μ = h0 .* p
    h1 ~ MvNormal(μ, σ)
end

m6_7 = sample(model_m6_7(d.h0, d.treatment, d.fungus, d.h1), NUTS(), 1000)
m6_7_df = DataFrame(m6_7)
precis(m6_7_df)
# -

# Code 6.17

# +
@model function model_m6_8(h0, treatment, h1)
    a ~ LogNormal(0, 0.2)
    bt ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    p = @. a + bt*treatment
    μ = h0 .* p
    h1 ~ MvNormal(μ, σ)
end

m6_8 = sample(model_m6_8(d.h0, d.treatment, d.h1), NUTS(), 1000)
m6_8_df = DataFrame(m6_8)
precis(m6_8_df)
# -

# Code 6.18

plant_dag = Dagitty.DAG(:H₀ => :H₁, :F => :H₁, :T => :F)
drawdag(plant_dag, [2, 0, 1, 3], [0, 0, 0, 0])

# Code 6.19

implied_conditional_independencies_min(plant_dag)

# Code 6.20

# +
Random.seed!(70)
# number of plants
N = 1000
h0 = rand(Normal(10, 2), N)
treatment = repeat(0:1, inner=div(N, 2))
M = rand(Bernoulli(), N)
fungus = [
    rand(Binomial(1, 0.5 - treat*0.4 + 0.4 * m)) 
    for (treat, m) ∈ zip(treatment, M)
]
h1 = h0 .+ rand(MvNormal(5 .+ 3 .* M, 1))

d2 = DataFrame(:h0 => h0, :h1 => h1, :treatment => treatment, :fungus => fungus)
precis(d2)
# -

m6_7 = sample(model_m6_7(d2.h0, d2.treatment, d2.fungus, d2.h1), NUTS(), 1000)
precis(DataFrame(m6_7))

m6_8 = sample(model_m6_8(d2.h0, d2.treatment, d2.h1), NUTS(), 1000)
precis(DataFrame(m6_8))

# # 6.3 Collider bias

# Code 6.21

d = sim_happiness(seed=1977, n_years=1000)
precis(d)

# +
d_m = d[d.married .== 1,[:age,:happiness]]
d_u = d[d.married .== 0,[:age,:happiness]]

scatter(d_m.age, d_m.happiness; label="married", xlab="age", ylab="happiness")
scatter!(d_u.age, d_u.happiness; c=:white)
# -

# Code 6.22

d2 = d[d.age .> 17,:]
d2[!,:A] = @. (d2.age - 18) / (65-18);

# Code 6.23

d2[!,:mid] = d2.married .+ 1;

# +
@model function model_m6_9(mid, A, happiness)
    a ~ MvNormal([0, 0], 1)
    bA ~ Normal(0, 2)
    μ = a[mid] .+ bA .* A
    σ ~ Exponential(1)
    happiness ~ MvNormal(μ, σ)
end

m6_9 = sample(model_m6_9(d2.mid, d2.A, d2.happiness), NUTS(), 1000)
m6_9_df = DataFrame(m6_9)
precis(m6_9_df)
# -

# Code 6.24

# +
@model function model_m6_10(A, happiness)
    a ~ Normal()
    bA ~ Normal(0, 2)
    μ = a .+ bA .* A
    σ ~ Exponential(1)
    happiness ~ MvNormal(μ, σ)
end

m6_10 = sample(model_m6_10(d2.A, d2.happiness), NUTS(), 1000)
m6_10_df = DataFrame(m6_10)
precis(m6_10_df)
# -

# Code 6.25

N = 200
b_GP = 1
b_GC = 0
b_PC = 1
b_U = 2;

# Code 6.26

Random.seed!(6)
U = 2 .* rand(Bernoulli(), N) .- 1
G = rand(Normal(), N)
P = rand(MvNormal(@. b_GP*G + b_U*U))
C = rand(MvNormal(@. b_PC*P + b_GC*G + b_U*U))
d = DataFrame(:C => C, :P => P, :G => G, :U => U);

# Code 6.27

# +
@model function model_m6_11(P, G, C)
    a ~ Normal()
    b_PC ~ Normal()
    b_GC ~ Normal()
    μ = @. a + b_PC*P + b_GC*G
    σ ~ Exponential(1)
    C ~ MvNormal(μ, σ)
end

m6_11 = sample(model_m6_11(d.P, d.G, d.C), NUTS(), 1000)
m6_11_df = DataFrame(m6_11)
precis(m6_11_df)
# -

# Code 6.28

# +
@model function model_m6_12(P, G, U, C)
    a ~ Normal()
    b_PC ~ Normal()
    b_GC ~ Normal()
    b_U ~ Normal()
    μ = @. a + b_PC*P + b_GC*G + b_U*U
    σ ~ Exponential(1)
    C ~ MvNormal(μ, σ)
end

m6_12 = sample(model_m6_12(d.P, d.G, d.U, d.C), NUTS(), 1000)
m6_12_df = DataFrame(m6_12)
precis(m6_12_df)
# -

# # 6.4 Confronting confounding

# Code 6.29

# +
dag_61 = Dagitty.DAG(
    :X => :Y,
    :U => :X, :A => :U,
    :A => :C, :C => :Y,
    :U => :B, :C => :B,
)

all_backdoor_adjustment_sets(dag_61, :X, :Y)
# -

# Code 6.30

dag_62 = Dagitty.DAG(
    :A => :D,
    :A => :M, :M => :D,
    :S => :A, :S => :M,
    :S => :W, :W => :D,
)
all_backdoor_adjustment_sets(dag_62, :W, :D)

implied_conditional_independencies_min(dag_62)


