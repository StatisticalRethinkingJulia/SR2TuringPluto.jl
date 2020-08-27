cd(@__DIR__)
using DrWatson
@quickactivate "StatReth"

# %%
using DataFrames
using CSV
using Distributions
using Turing
using Plots
# using Statistics

include(srcdir("quap.jl"))
include(srcdir("tools.jl"))

# %% 4.64
d = DataFrame(CSV.File(datadir("exp_raw/Howell_1.csv")))

# %% 4.65, 6.66
d.weight_s = (d.weight .- mean(d.weight)) / std(d.weight)

f_parabola(weight_s, a, b1, b2) = a + b1 * weight_s + b2 * weight_s^2

@model function parabola(weight_s, heights)
    a ~ Normal(178, 20)
    b1 ~ LogNormal(0, 1)
    b2 ~ Normal(0, 1)
    σ ~ Uniform(0, 50)
    μ = f_parabola.(weight_s, a, b1, b2)
    heights ~ MvNormal(μ, σ)
end

q4_5 = quap(parabola(d.weight_s, d.height), NelderMead())

# precis(m4_5)

# %% 4.67
weight_seq = range(-2.2, 2, length = 30)
post = DataFrame(rand(q4_5.distr, 1_000)', q4_5.params)
mu = f_parabola.(weight_seq, post.a', post.b1', post.b2')
sim = rand.(Normal.(mu, post.σ'))

mu = meanlowerupper(mu)
sim = meanlowerupper(sim)

# %% 4.68
scatter(d.weight_s, d.height, ms = 3, alpha = 0.7, legend = false)
plot!(weight_seq, mu.mean, ribbon = (mu.mean .- mu.lower, mu.upper .- mu.mean))
plot!(weight_seq, sim.lower, fillrange = sim.upper, alpha = 0.3, linealpha = 0.0, c = 2)

# %% 4.69
f_cube(weight_s, a, b1, b2, b3) = a + b1 * weight_s + b2 * weight_s^2 + b3 * weight_s^3

@model function cube(weight_s, heights)
    a ~ Normal(178, 20)
    b1 ~ LogNormal(0, 1)
    b2 ~ Normal(0, 10)
    b3 ~ Normal(0, 10)
    σ ~ Uniform(0, 50)
    μ = f_cube.(weight_s, a, b1, b2, b3)
    heights ~ MvNormal(μ, σ)
end

q4_6 = quap(cube(d.weight_s, d.height), NelderMead())

# %% 4.70, 4.71
weight_seq = range(-2.2, 2, length = 30)
post = DataFrame(rand(q4_6.distr, 1_000)', q4_6.params)
mu = f_cube.(weight_seq, post.a', post.b1', post.b2', post.b3') |> meanlowerupper
sim = rand.(Normal.(mu.raw, post.σ')) |> meanlowerupper

weight_seq_rescaled = weight_seq .* std(d.weight) .+ mean(d.weight)
scatter(d.weight, d.height, ms = 3, alpha = 0.7, legend = false)
plot!(weight_seq_rescaled, mu.mean, ribbon = (mu.mean .- mu.lower, mu.upper .- mu.mean))
plot!(weight_seq_rescaled, sim.lower, fillrange = sim.upper, alpha = 0.3, la = 0.0, c = 2)

# %% 4.72
d = DataFrame!(CSV.File(datadir("exp_raw/cherry_blossoms.csv"), missingstring = "NA"))

# precis(d)

scatter(d.year, d.doy)

# %% 4.73
d2 = dropmissing(d, :doy)

num_knots = 15
knot_list = quantile(d2.year, range(0, 1, length = num_knots))

# %% 4.74, 4.75
using BSplines: BSplineBasis, basismatrix

basis = BSplineBasis(4, knot_list)
B = basismatrix(basis, d2.year)

plot(legend = false, xlabel = "year", ylabel = "basis value")
for y in eachcol(B)
    plot!(d2.year, y)
end
plot!()

# %% 4.76
@model function spline(D, B = B)
    α ~ Normal(100, 10)
    w ~ filldist(Normal(0, 10), size(B, 2))
    σ ~ Exponential(1)
    μ = α .+ B * w
    D ~ MvNormal(μ, σ)
    return μ
end

q4_7 = quap(spline(d2.doy))

# %% 4.77
w_str = ["w[$i]" for i in 1:length(basis)]
post = DataFrame(rand(q4_7.distr, 1000)', ["α"; w_str; "σ"])

w = mean.(eachcol(post[:, w_str]))              # either
w = [mean(post[:, col]) for col in w_str]       # or

plot(legend = false, xlabel = "year", ylabel = "basis * weight")
for y in eachcol(B .* w')
    plot!(d2.year, y)
end
plot!()

# %% 4.78
mu = post.α' .+ B * Array(post[!, w_str])'
mu = meanlowerupper(mu)

scatter(d2.year, d2.doy, alpha = 0.3)
plot!(d2.year, mu.mean, ribbon = (mu.mean .- mu.lower, mu.upper .- mu.mean))

# %% 4.79
@model function spline(D, B = B)
    α ~ Normal(100, 10)
    w ~ filldist(Normal(0, 10), size(B, 2))
    σ ~ Exponential(1)
    μ = [α + sum(Brow .* w) for Brow in eachrow(B)]
    D ~ MvNormal(μ, σ)
end

q4_7alt = quap(spline(d2.doy))
