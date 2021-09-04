using MCMCChains
using Distributions
using ParetoSmooth

samples = randn(100, 1, 1)
data = randn(50)
chain = Chains(samples)

compute_loglike(μ, data) = logpdf(Normal(μ, 1), data)
compute_loglike(μ, σ, data) = logpdf(Normal(μ, σ), data)

pll1 = pointwise_log_likelihoods(compute_loglike, chain, data)
size(pll1) |> display
