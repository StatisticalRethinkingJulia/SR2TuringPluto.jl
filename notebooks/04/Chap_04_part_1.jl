using Distributions
using StatsPlots
using LaTeXStrings
using CSV
using DataFrames
using StatisticalRethinking
using LinearAlgebra
using Logging
using Random
using Turing

# setting default attributes for plots
default(label=false)
Logging.disable_logging(Logging.Warn);
# -

# # 4.1 Why normal distributions are normal

# Code 4.1

n = rand(Uniform(-1, 1), 1000, 16);
pos = sum.(eachrow(n));
density(pos)

# Code 4.2

prod(1 .+ rand(Uniform(0, .1), 12))

# Code 4.3

# +
u = Uniform(0, .1)
growth = prod.(eachrow(1 .+ rand(u, 10000, 12)));

density(growth; label="density")
# overlay normal distribution
μ = mean(growth)
σ = std(growth)
plot!(Normal(μ, σ); label="Normal")
# -

# Code 4.4

big = prod.(eachrow(1 .+ rand(Uniform(0, 0.5), 10000, 12)));
small = prod.(eachrow(1 .+ rand(Uniform(0, 0.01), 10000, 12)));
density([big, small]; layout=(2, 1))

# Code 4.5

density(log.(big))

# # 4.2 A language for describing models

# Code 4.6

w = 6
n = 9
p_grid = range(0, 1; length=100)
bin_dens = [pdf(Binomial(n, p), w) for p in p_grid]
uni_dens = [pdf(Uniform(0, 1), p) for p in p_grid];
posterior = bin_dens .* uni_dens
posterior /= sum(posterior);

# # 4.3 Gaussian model of height

# Code 4.7

d = DataFrame(CSV.File("data/Howell1.csv"));

# Code 4.8

describe(d)

# Code 4.9

precis(d)

# Code 4.10

d.height

# Code 4.11

d2 = d[d.age .>= 18,:];

# Code 4.12

plot(Normal(178, 20); xlim=(100, 250))

# Code 4.13

plot(Uniform(0, 50), xlim=(-10, 60), ylim=(0, 0.03))

# Code 4.14

# +
size = 10_000
sample_μ = rand(Normal(178, 20), size)
sample_σ = rand(Uniform(0, 50), size);
prior_h = [rand(Normal(μ, σ)) for (μ, σ) in zip(sample_μ, sample_σ)];

p1 = density(sample_μ; title="μ")
p2 = density(sample_σ; title="σ")
p3 = density(prior_h; title="prior_h")

plot(p1, p2, p3, layout=@layout [a b; c])
# -

# Code 4.15

# +
sample_μ = rand(Normal(178, 100), size)
prior_h = [rand(Normal(μ, σ)) for (μ, σ) in zip(sample_μ, sample_σ)];

density(prior_h)
vline!([0, 272])
# -

# Code 4.16

# +
μ_list = range(150, 160; length=100)
σ_list = range(7, 9; length=100)

log_likelihood = [
    sum(logpdf(Normal(μ, σ), d2.height))
    for μ ∈ μ_list, σ ∈ σ_list
]
log_prod = log_likelihood .+ [
    logpdf(Normal(178, 20), μ) + logpdf(Uniform(0, 50), σ)
    for μ ∈ μ_list, σ ∈ σ_list
];

max_prod = maximum(log_prod)
prob = @. exp(log_prod - max_prod);
# -

# Code 4.17

# note the transposition, that's due to Julia matrix order
contour(μ_list, σ_list, prob')

# Code 4.18

heatmap(μ_list, σ_list, prob')

# Code 4.19

indices = collect(Iterators.product(1:length(μ_list), 1:length(σ_list)));
sample_idx = wsample(vec(indices), vec(prob), 10_000; replace=true)
sample_μ = μ_list[first.(sample_idx)]
sample_σ = σ_list[last.(sample_idx)];

# Code 4.20

scatter(sample_μ, sample_σ; alpha=0.1)

# Code 4.21

p1 = density(sample_μ)
p2 = density(sample_σ)
plot(p1, p2, layout=(2,1))

# Code 4.22

println(hpdi(sample_μ, alpha=0.11))
println(hpdi(sample_σ, alpha=0.11))

# Code 4.23

d3 = sample(d2.height, 20);

# Code 4.24

# +
μ_list = range(150, 170; length=200)
σ_list = range(4, 20; length=200)

log_likelihood = [
    sum(logpdf(Normal(μ, σ), d3))
    for μ ∈ μ_list, σ ∈ σ_list
]
log_prod = log_likelihood .+ [
    logpdf(Normal(178, 20), μ) + logpdf(Uniform(0, 50), σ)
    for μ ∈ μ_list, σ ∈ σ_list
]

max_prod = maximum(log_prod)
prob2 = @. exp(log_prod - max_prod)

indices = collect(Iterators.product(1:length(μ_list), 1:length(σ_list)));
sample2_idx = wsample(vec(indices), vec(prob2), 10_000; replace=true)
sample2_μ = μ_list[first.(sample2_idx)]
sample2_σ = σ_list[last.(sample2_idx)]

scatter(sample2_μ, sample2_σ; alpha=0.1)
# -

# Code 4.25

density(sample2_σ)
μ = mean(sample2_σ)
σ = std(sample2_σ)
plot!(Normal(μ, σ); label="Normal")

# Code 4.26

d = DataFrame(CSV.File("data/Howell1.csv"));
d2 = d[d.age .>= 18,:];

# Code 4.27

@model function model_height(height)
    μ ~ Normal(178, 20)
    σ ~ Uniform(0, 50)
    height ~ Normal(μ, σ)
end

# Code 4.28

m4_1 = sample(model_height(d2.height), NUTS(), 1000)

# Code 4.29

display.(describe(m4_1; q=[0.055, 0.945]));

# Code 4.30

init_vals = [mean(d2.height), std(d2.height)]
chain = sample(model_height(d2.height), NUTS(), 1000, init_theta=init_vals)

# Code 4.31

# +
@model function model_height(height)
    μ ~ Normal(178, 0.1)
    σ ~ Uniform(0, 50)
    height ~ Normal(μ, σ)
end

m4_2 = sample(model_height(d2.height), NUTS(), 1000)
display.(describe(m4_2; q=[0.055, 0.945]));
# -

# Code 4.32

cov(hcat(m4_1[:μ], m4_1[:σ]))

# Code 4.33

c = cov(hcat(m4_1[:μ], m4_1[:σ]))
cov2cor(c, diag(c))

# Code 4.34

# resetrange is needed due to bug in MCMCChains: https://github.com/TuringLang/MCMCChains.jl/issues/324
# once it will be fixed, direct sampling from the chain will be possible
samp_chain = sample(resetrange(m4_1), 10_000)
samp_df = DataFrame(samp_chain)
first(samp_df, 5)

# Code 4.35

precis(samp_df)

# Code 4.36

data = hcat(m4_1[:μ], m4_1[:σ])
μ = mean(data, dims=1)
σ = cov(data)
mvn = MvNormal(vec(μ), σ)
post = rand(mvn, 10_000);
print(mvn)

# # 4.4 Linear predictions

# Code 4.37

d = DataFrame(CSV.File("data/Howell1.csv"));
d2 = d[d.age .>= 18,:];
# fancy way of doing scatter(d2.weight, d2.height)
@df d2 scatter(:weight, :height)

# Code 4.38

Random.seed!(2971)
N = 100
a = rand(Normal(178, 20), N)
b = rand(Normal(0, 10), N);

# Code 4.39

# +
p = hline([0, 272]; ylims=(-100, 400), xlabel="weight", ylabel="hegiht")
title!(L"\beta \sim \mathcal{N}(\mu=0,\sigma=10)")

x_mean = mean(d2.weight)
xlims = extrema(d2.weight)  # getting min and max in one pass

for (α, β) ∈ zip(a, b)
    plot!(x -> α + β * (x - x_mean); xlims=xlims, c=:black, alpha=0.3)
end
display(p)
# -

# Code 4.40

b = rand(LogNormal(0, 1), 10_000)
density(b, xlims=(0, 5), bandwidth=0.1)

# Code 4.41

Random.seed!(2971)
N = 100
a = rand(Normal(178, 20), N)
b = rand(LogNormal(0, 1), N);

# Code 4.42

# +
d = DataFrame(CSV.File("data/Howell1.csv"));
d2 = d[d.age .>= 18,:]
xbar = mean(d2.weight)

@model function height_regr_model(weight, height)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    μ = @. a + b * (weight - xbar)
    σ ~ Uniform(0, 50)
    height ~ MvNormal(μ, σ)
end

m4_3 = sample(height_regr_model(d2.weight, d2.height), NUTS(), 1000)
m4_3 = resetrange(m4_3);
# -

# Code 4.43

# +
@model function height_regr_model_exp(weight, height)
    a ~ Normal(178, 20)
    log_b ~ Normal(0, 1)
    μ = @. a + exp(log_b) * (weight - xbar)
    σ ~ Uniform(0, 50)
    height ~ MvNormal(μ, σ)
end

m4_3b = sample(height_regr_model_exp(d2.weight, d2.height), NUTS(), 1000);
# -

# Code 4.44

m4_3_df = DataFrame(m4_3)
precis(m4_3_df)

# Code 4.45

round.(cov(Matrix(m4_3_df)), digits=3)

# Code 4.46

# +
p = @df d2 scatter(:weight, :height; alpha=0.3)

chain = resetrange(m4_3)
samples = sample(chain, 1000)

a_map = mean(samples[:a])
b_map = mean(samples[:b])
plot!(x -> a_map + b_map*(x-xbar))
# -

# Code 4.47

post = sample(m4_3, 1000)
post_df = DataFrame(post)
post_df[1:5,:]

# Code 4.48

# +
N = 10
dN = d2[1:N,:]

@model function height_regr_model_N(weight, height)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    m_weight = mean(weight)
    μ = @. a + b * (weight - m_weight)
    σ ~ Uniform(0, 50)
    height ~ MvNormal(μ, σ)
end

mN = sample(height_regr_model(dN.weight, dN.height), NUTS(), 1000)
mN = resetrange(mN);
# -

# Code 4.49

# +
post = sample(mN, 20)
post_df = DataFrame(post);

xlims = extrema(d2.weight)
ylims = extrema(d2.height)
p = @df dN scatter(:weight, :height; xlims=xlims, ylims=ylims)
title!("N = $N"; xlab="weight", ylab="height")

x_mean = mean(dN.weight)
for (a, b) ∈ zip(post_df.a, post_df.b)
    plot!(x -> a + b * (x-x_mean); c="black", alpha=0.3)
end
display(p)
# -

# Code 4.50

post = sample(m4_3, 1000)
post_df = DataFrame(post)
μ_at_50 = @. post_df.a + post_df.b * (50 - xbar);

# Code 4.51

density(μ_at_50; lw=2, xlab=L"\mu|weight=50")

# Code 4.52

PI(μ_at_50; prob=0.89)

# Code 4.53

μ = StatisticalRethinking.link(post_df, [:a :b], d2.weight, xbar);
μ = hcat(μ...);
Base.size(μ), μ[1:5,1]

# Code 4.54

weight_seq = 25:70
μ = StatisticalRethinking.link(post_df, [:a :b], weight_seq, xbar);
μ = hcat(μ...);
Base.size(μ), μ[1:5,1]

# Code 4.55

p = plot()
for i in 1:100
    scatter!(weight_seq, μ[i,:]; c=:blue, alpha=0.2)
end
display(p)

# Code 4.56

μ_mean = mean.(eachcol(μ))
μ_PI = PI.(eachcol(μ))
μ_PI = vcat(μ_PI'...);

# Code 4.57

@df d2 scatter(:weight, :height; alpha=0.2, xlab="weight", ylab="height")
plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.3)

# Code 4.58

# +
post = sample(m4_3, 1000)
post = DataFrame(post)

weight_seq = 25:70
μ = map(w -> post.a + post.b * (w - xbar), weight_seq)
μ = hcat(μ...)
μ_mean = mean.(eachcol(μ))
μ_CI = PI.(eachcol(μ));
# -

# Code 4.59

sim_height = simulate(post, [:a, :b, :σ], weight_seq .- xbar);
Base.size(sim_height), sim_height[1:5,1]

# Code 4.60

height_PI = PI.(eachcol(sim_height))
height_PI = vcat(height_PI'...);

# Code 4.61

@df d2 scatter(:weight, :height; alpha=0.2, xlab="weight", ylab="height")
plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.3)
plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=height_PI, fillalpha=0.3)

# Code 4.62

post = sample(m4_3, 10_000)
post = DataFrame(post)
sim_height = simulate(post, [:a, :b, :σ], weight_seq .- xbar)
height_PI = PI.(eachcol(sim_height))
height_PI = vcat(height_PI'...);

@df d2 scatter(:weight, :height; alpha=0.2, xlab="weight", ylab="height")
plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.3)
plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=height_PI, fillalpha=0.3)

# Code 4.63

# +
post = sample(m4_3, 1000)
post = DataFrame(post)

sim_height = [
    [
        rand(Normal(a + b * (w - xbar), σ))
        for (a, b, σ) ∈ zip(post.a, post.b, post.σ)
    ]
    for w ∈ weight_seq
]
sim_height = hcat(sim_height...)

height_PI = PI.(eachcol(sim_height));
height_PI = vcat(height_PI'...);
# -

@df d2 scatter(:weight, :height; alpha=0.2, xlab="weight", ylab="height")
plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.3)
plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=height_PI, fillalpha=0.3)

# # 4.5 Curves from lines

# Code 4.64

d = DataFrame(CSV.File("data/Howell1.csv"))
scatter(d.weight, d.height; alpha=0.3)

# Code 4.65

# +
d[!, :weight_s] = standardize(ZScoreTransform, d.weight)
d[!, :weight_s2] = d.weight_s.^2;

@model function height_regr_m2(weight_s, weight_s2, height)
    a ~ Normal(178, 20)
    b1 ~ LogNormal(0, 1)
    b2 ~ Normal(0, 1)
    μ = @. a + b1 * weight_s + b2 * weight_s2
    σ ~ Uniform(0, 50)
    height ~ MvNormal(μ, σ)
end

m4_5 = sample(height_regr_m2(d.weight_s, d.weight_s2, d.height), NUTS(), 1000)
m4_5 = resetrange(m4_5)
m4_5_df = DataFrame(m4_5);
# -

# Code 4.66

precis(m4_5_df)

# Code 4.67

# +
Random.seed!(1)
df = sample(m4_5_df, 1000)
weight_seq = range(-2.2, 2; length=30)

# explicit logic of link
mu = [
    df.a + df.b1 * w_s + df.b2 * w_s^2
    for w_s ∈ weight_seq
]

mu = hcat(mu...)
mu_mean = mean.(eachcol(mu))
mu_PI = PI.(eachcol(mu))
mu_PI = vcat(mu_PI'...)

# explicit logic of sim
sim_height = [
    [
        rand(Normal(row.a + row.b1 * w_s + row.b2 * w_s^2, row.σ))
        for row ∈ eachrow(df)
    ]
    for w_s ∈ weight_seq
]
sim_height = hcat(sim_height...);

height_PI = PI.(eachcol(sim_height))
height_PI = vcat(height_PI'...);
# -

# Code 4.68

p_square = @df d scatter(:weight_s, :height; alpha=0.3, title="Square poly")
plot!(weight_seq, mu_mean; c=:black)
plot!(weight_seq, [mu_mean mu_mean]; c=:black, fillrange=mu_PI, fillalpha=0.3)
plot!(weight_seq, [mu_mean mu_mean]; c=:black, fillrange=height_PI, fillalpha=0.3)

# Code 4.69

# +
d[!, :weight_s3] = d.weight_s.^3;

@model function height_regr_m3(weight_s, weight_s2, weight_s3, height)
    a ~ Normal(178, 20)
    b1 ~ LogNormal(0, 1)
    b2 ~ Normal(0, 1)
    b3 ~ Normal(0, 1)
    μ = @. a + b1 * weight_s + b2 * weight_s2 + b3 * weight_s3
    σ ~ Uniform(0, 50)
    height ~ MvNormal(μ, σ)
end

m4_6 = sample(height_regr_m3(d.weight_s, d.weight_s2, d.weight_s3, d.height), NUTS(), 1000)
m4_6 = resetrange(m4_6)
m4_6_df = DataFrame(m4_6)
precis(m4_6_df)

# +
Random.seed!(1)
df = sample(m4_6_df, 1000)
weight_seq = range(-2.2, 2; length=30)

# explicit logic of link
mu = [
    df.a + df.b1 * w_s + df.b2 * w_s^2 + df.b3 * w_s^3
    for w_s ∈ weight_seq
]

mu = hcat(mu...)
mu_mean = mean.(eachcol(mu))
mu_PI = PI.(eachcol(mu))
mu_PI = vcat(mu_PI'...)

# explicit logic of sim
sim_height = [
    [
        rand(Normal(row.a + row.b1 * w_s + row.b2 * w_s^2 + row.b3 * w_s^3, row.σ))
        for row ∈ eachrow(df)
    ]
    for w_s ∈ weight_seq
]
sim_height = hcat(sim_height...);

height_PI = PI.(eachcol(sim_height))
height_PI = vcat(height_PI'...);

# +
p_cubic = @df d scatter(:weight_s, :height; alpha=0.3, title="Cubic poly")
plot!(weight_seq, mu_mean; c=:black)
plot!(weight_seq, [mu_mean mu_mean]; c=:black, fillrange=mu_PI, fillalpha=0.3)
plot!(weight_seq, [mu_mean mu_mean]; c=:black, fillrange=height_PI, fillalpha=0.3)

plot(p_square, p_cubic; layout=(1, 2))
# -

# Code 4.70 and 4.71
#
# Looks like Julia plots don't support change of ticks proposed in the book.
# But much more natural way will be to remap values we're plotting back to the original scale.
# Example of this is below.

# +
μ = mean(d.weight)
σ = std(d.weight)
w = @. d.weight_s * σ + μ
scatter(w, d.height; alpha=0.3)

w_s = @. weight_seq * σ + μ
plot!(w_s, mu_mean; c=:black)
plot!(w_s, [mu_mean mu_mean]; c=:black, fillrange=mu_PI, fillalpha=0.3)
plot!(w_s, [mu_mean mu_mean]; c=:black, fillrange=height_PI, fillalpha=0.3)
# -

# Code 4.72

d = DataFrame(CSV.File("data/cherry_blossoms.csv", missingstring="NA"))
precis(d)

# Code 4.73

d2 = d[completecases(d[!,[:doy]]),:]
d2 = disallowmissing(d2[!,[:year,:doy]])
num_knots = 15
knots_list = quantile(d2.year, range(0, 1; length=num_knots));

# Code 4.74

# +
using BSplines

basis = BSplineBasis(3, knots_list)
# -

# Code 4.75

p1 = plot(basis)
scatter!(knots_list, repeat([1], num_knots); xlab="year", ylab="basis", legend=false)

# Code 4.76
#
# This way of calucalting bsplines is slightly slower, than shown in the book (with pre-calculated matrix) but it is much cleaner in my perspective.
#
# You can do comparison yourself by precalculating spline matrix outside of the model and do matrix multiplication in the model instead of spline evialutaion. Example of doing this is at code block 4.79

# +
@model function model_splines(year, doy)
    w ~ MvNormal(zeros(length(basis)), 1)
    a ~ Normal(100, 10)
    s = Spline(basis, w)
    μ = a .+ s.(year)
    σ ~ Exponential(1)
    doy ~ MvNormal(μ, σ)
end

m4_7 = sample(model_splines(d2.year, d2.doy), NUTS(0.65; init_ϵ = 9.765625e-5), 1000)
# -

# Code 4.77

# +
post = DataFrame(m4_7)

# convert columns w[*] into single column w
w_df = DataFrames.select(post, r"w")
post = DataFrames.select(post, Not(r"w"))
post[!,:w] = Vector.(eachrow(w_df))

# vector of 16 average w values
w_mean = mean.(eachcol(w_df))
p2 = plot(basis .* w_mean)
scatter!(knots_list, repeat([1], num_knots); xlab="year", ylab="basis × weight")
# -

# Code 4.78

# +
# explicit link logic
μ = [
    row.a .+ Spline(basis, row.w).(d2.year)
    for row ∈ eachrow(post)
]
μ = hcat(μ...);

μ_PI = PI.(eachrow(μ))
μ_PI = vcat(μ_PI'...);

p3 = @df d2 scatter(:year, :doy; alpha=0.3)
μ_mean = mean.(eachrow(μ_PI))
plot!(d2.year, [μ_mean, μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.3, alpha=0)
# -

plot(p1, p2, p3; layout=(3, 1))

# Code 4.79
#
# How to build the model with explicit spline matrix calculation

# +
basis = BSplineBasis(3, knots_list)

# list of splines with 1 only at corresponding basis index
splines = [
    Spline(basis, [float(idx == knot) for idx ∈ 1:length(basis)])
    for knot ∈ 1:length(basis)
]

# calculate each spline for every year. Resulting matrix B is 827x16
B = [
    map(s -> s(year), splines)
    for year in d2.year
]
B = vcat(B'...);


# do not need years parameter anymore, all the information is in B matrix
@model function model_splines_matrix(doy)
    w ~ MvNormal(zeros(length(basis)), 1)
    a ~ Normal(100, 10)
    μ = a .+ B * w
    σ ~ Exponential(1)
    doy ~ MvNormal(μ, σ)
end

m4_7alt = sample(model_splines_matrix(d2.doy), NUTS(0.65; init_ϵ = 0.0001953125), 1000)
# -


