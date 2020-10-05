# Clip-05-01-27.jl

using DrWatson
@quickactivate "StatReth"
using Turing

include(srcdir("quap.jl"))

# %% 5.1, 5.2
d = DataFrame(CSV.File(datadir("exp_raw/WaffleDivorce.csv")))
d.D = zscore(d.Divorce)
d.M = zscore(d.Marriage)
d.A = zscore(d.MedianAgeMarriage)

std(d.MedianAgeMarriage)

# %% 5.3
@model function divorce_A(A, D)
    a ~ Normal(0, 0.2)
    bA ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bA)
    D ~ MvNormal(μ, σ)
end

m5_1 = divorce_A(d.A, d.D)
prior = sample(m5_1, Prior(), 50) |> DataFrame

# %% 5.4
x = -2:0.1:2
plot()
for r in eachrow(prior)
    p = lin(r.a, x, r.bA)
    plot!(x, p, color = :black, alpha = 0.4)
end
plot!(legend = false)

# %% 5.5
q5_1 = quap(m5_1)
post = DataFrame(rand(q5_1.distr, 1000)', q5_1.params)

A_seq = range(-3, 3.2, length = 30)
mu = lin(post.a', A_seq, post.bA') |> meanlowerupper

scatter(d.A, d.D, alpha = 0.4, legend = false)
plot!(A_seq, mu.mean, ribbon = (mu.mean .- mu.lower, mu.upper .- mu.mean))

# %% 5.6
@model function divorce_M(M, D)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, M, bM)
    D ~ MvNormal(μ, σ)
end

q5_2 = quap(divorce_M(d.M, d.D))
post = DataFrame(rand(q5_2.distr, 1000)', q5_2.params)

M_seq = range(-3, 3.2, length = 30)
mu = lin(post.a', M_seq, post.bM') |> meanlowerupper

scatter(d.M, d.D, alpha = 0.4, legend = false)
plot!(M_seq, mu.mean, ribbon = (mu.mean .- mu.lower, mu.upper .- mu.mean))

# %% 5.7 - 5.9
missing  # TODO

# %% 5.10
@model function divorce_A_M(A, M, D)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    bA ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bA, M, bM)
    D ~ MvNormal(μ, σ)
end

m5_3 = divorce_A_M(d.A, d.M, d.D)
q5_3 = quap(m5_3)

# %% 5.11
# Well, this ain't pretty. But before I have something reasonable, I at least wanted to
# something, even if it's messily thrown together. TODO

x = [
    q5_1.coef.bA,       # bA
    NaN,
    q5_3.coef.bA,
    NaN,                # Spaceholder
    NaN,                # bM
    q5_2.coef.bM,
    q5_3.coef.bM,
]
xerr = sqrt.([
    q5_1.vcov[2, 2],    # bA
    NaN,
    q5_3.vcov[3, 3],
    NaN,                # Spaceholder
    NaN,                # bM
    q5_2.vcov[2, 2],
    q5_3.vcov[2, 2],
])
ylab = [
    "bA m5.1",
    "bA m5.2",
    "bA m5.3",
    "",
    "bM m5.1",
    "bM m5.2",
    "bM m5.3",
]

scatter(x, ylab, xerr = xerr, legend = false)
vline!([0])

# %% 5.12
N = 50
age = randn(N)
mar = rand.(Normal.(-age))
div = rand.(Normal.(age))

# %% 5.13, 5.14
@model function divorce_AM(A, M)
    a ~ Normal(0, 0.2)
    bAM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bAM)
    M ~ MvNormal(μ, σ)
end

sort!(d, :A)  # for nice looking graphs

q5_4 = quap(divorce_AM(d.A, d.M))
post = DataFrame(rand(q5_4.distr, 1000)', q5_4.params)
mu = lin(post.a', d.A, post.bAM') |> meanlowerupper
resid = d.M .- mu.mean

scatter(d.A, d.M, legend = false)
plot!(d.A, mu.mean, ribbon = (mu.mean .- mu.lower, mu.upper .- mu.mean))

# %% 5.15, 5.16
post = DataFrame(rand(q5_3.distr, 1000)', q5_3.params)
μ = lin(post.a', d.A, post.bA', d.M, post.bM') |> meanlowerupper

scatter(d.D, μ.mean, yerror = (μ.mean .- μ.lower, μ.upper .- μ.mean), legend = false)
plot!(identity, range(extrema(d.D)..., length = 10))
plot!(xlabel = "Observed divorce", ylabel = "Predicted divorce")

# %% 5.17
xy = [d.D μ.mean]
for loc in ("ID", "UT", "ME", "RI")
    coord = xy[d.Loc .== loc, :]
    annotate!(coord..., loc, :bottom)
end
plot!()

# %% 5.18
N = 100
x_real = randn(N)
x_spur = rand.(Normal.(x_real))
y = rand.(Normal.(x_real))
d = DataFrame((; y, x_real, x_spur))  # the semicolon turns it into a NamedTuple which then
                                      # gives the DataFrame the names of the columns

@df d corrplot([:y :x_real :x_spur])  # either
corrplot(Matrix(d), label = names(d)) # or

# %% 5.19
d = DataFrame(CSV.File(datadir("exp_raw/WaffleDivorce.csv")))
d.D = zscore(d.Divorce)
d.M = zscore(d.Marriage)
z_age, unz_age = zscore_transform(d.MedianAgeMarriage)
d.A = z_age(d.MedianAgeMarriage)

@model function divorce_19(A, M, D)
    # A -> M
    aM ~ Normal(0, 0.2)
    bAM ~ Normal(0, 0.5)
    σ_M ~ Exponential(1)
    μ_M = lin(aM, A, bAM)
    M ~ MvNormal(μ_M, σ_M)
    # A -> D <- M
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    bA ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bA, M, bM)
    D ~ MvNormal(μ, σ)
end

q5_3_A = quap(divorce_19(d.A, d.M, d.D))

# %% 5.20 - 5.22
A_seq = range(-2, 2, length = 30)
post = DataFrame(rand(q5_3_A.distr, 1000)', q5_3_A.params)

μM = lin(post.aM', A_seq, post.bAM')
M = rand.(Normal.(μM, post.σ_M'))
μD = lin(post.a', A_seq, post.bA', M, post.bM')
D = rand.(Normal.(μD, post.σ')) |> meanlowerupper

plot(A_seq, D.mean, ribbon = (D.mean .- D.lower, D.upper .- D.mean))
plot!(legend = false, xlabel = "manipulated A", ylabel = "counterfactual D")

# %% 5.23
# I use zscore_transform from the tools.jl file so I don't have to hardcode the mean and
# stddev of d.MedianAgeMarriage.
M = rand.(Normal.(lin(post.aM', (20.0, 30.0) |> z_age, post.bAM'), post.σ_M'))
μD = lin(post.a', (20.0, 30.0) |> z_age, post.bA', M, post.bM') |> meanlowerupper
μD.mean[2] - μD.mean[1]
# Eh, yes, 4.5 stddev are kinda big but so is a change in the median age of marriage from
# 20 to 30 years, clocking in at 8 stddevs. I mean, think about it, a median age of
# marriage of 20 years?!

# %% 5.24
M_seq = range(-2, 2, length = 30)
A_seq = zeros(length(M_seq))
M = rand.(Normal.(lin(post.aM', A_seq, post.bAM'), post.σ_M'))
μD = lin(post.a', A_seq, post.bA', M_seq, post.bM')
D = rand.(Normal.(μD, post.σ')) |> meanlowerupper

plot(M_seq, D.mean, ribbon = (D.mean .- D.lower, D.upper .- D.mean))
plot!(legend = false, xlabel = "manipulated M", ylabel = "counterfactual D")

# %% 5.25 - 5.27
# This is a little superfluous since we aren't hiding any details but here we go.
A_seq = range(-2, 2, length = 30)
post = DataFrame(rand(q5_3_A.distr, 1000)', q5_3_A.params)
M_sim = rand.(Normal.(post.aM' .+ A_seq .* post.bAM', post.σ_M'))
D_sim = rand.(Normal.(post.a' .+ A_seq .* post.bA', post.σ'))

# End of clip-05-01-27.jl
