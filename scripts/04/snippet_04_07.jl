using DrWatson
@quickactivate "StatReth"

# %%
using DataFrames
using CSV
using StatsBase
using Distributions
using StatsPlots

# %% 4.7 - 4.11
d = DataFrame(CSV.File(datadir("exp_raw/Howell_1.csv")))
# precis(d)
d.height

d2 = filter(row -> row.age >= 18, d)    # either
d2 = d[d.age .>= 18, :]                 # or

# %% 4.12, 4.13
plot(p -> pdf(Normal(178, 20), p), 100:250)
plot(p -> pdf(Uniform(0, 50), p), -10:60)

# %% 4.14
N = 10_000
d_μ = Normal(178, 20)
sample_mu = rand(d_μ, N)
d_σ = Uniform(0, 50)
sample_sigma = rand(d_σ, N)

prior_h = rand.(Normal.(sample_mu, sample_sigma))
density(prior_h)

# %% 4.15
d_μ_100 = Normal(178, 100)
sample_μ_100 = rand(d_μ_100, N)

sample_height_100 = rand.(Normal.(sample_μ_100, sample_sigma))
# or is
# mix = MixtureModel(Normal.(sample_μ_100, sample_sigma))
# sample_height_100 = rand(mix, N)
# more correct?
density(sample_height_100)

# %% 4.16, 4.17
μ_grid = range(150, 160, length = 100)
σ_grid = range(7, 9, length = 100)
post = [(; :μ => μ, :σ => σ) for μ in μ_grid, σ in σ_grid]  # this is `post`

ll = [loglikelihood(Normal(μ, σ), d2.height) for (μ, σ) in post]
prior = [logpdf(d_μ, μ) + logpdf(d_σ, σ) for (μ, σ) in post]
lposterior = ll .+ prior
lposterior .-= maximum(lposterior)
posterior = exp.(lposterior)

wireframe(μ_grid, σ_grid, posterior)

# %% 4.19 - 4.22
samples = sample(post, pweights(posterior), N) |> DataFrame

scatter(samples.μ, samples.σ, alpha = 0.05, ms = 5, xlabel = "μ", ylabel = "σ")

density(samples.μ)
density(samples.σ)

quantile(samples.μ, (0.1, 0.9))
quantile(samples.σ, (0.1, 0.9))

# %% 4.23 - 4.25
d3_height = d2.height[1:20]
μ_grid = range(140, 160, length = 200)
σ_grid = range(4, 20, length = 200)
post = [(μ = μ, σ = σ) for μ in μ_grid, σ in σ_grid]

ll = [loglikelihood(Normal(μ, σ), d3_height) for (μ, σ) in post]
prior = [logpdf(d_μ, μ) + logpdf(d_σ, σ) for (μ, σ) in post]
lposterior = ll .+ prior
lposterior .-= maximum(lposterior)
posterior = exp.(lposterior)

wireframe(μ_grid, σ_grid, posterior)

samples2 = sample(post, pweights(posterior), N) |> DataFrame

scatter(samples2.μ, samples2.σ, alpha = 0.05, ms = 5)
density(samples2.μ)
density(samples2.σ)
