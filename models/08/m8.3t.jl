# ### m8.3t

using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking
Turing.turnprogress(false)

# Turing model
@model ppl8_3(y) = begin
    α ~ Normal(1, 10)
    σ ~ truncated(Cauchy(0, 1), 0, Inf)

    y .~ Normal(α, σ)
end

y = [-1,1]

# Sample

m8_3t = ppl8_3(y)
nchains = 4; sampler = NUTS(0.65); nsamples=2000
chns8_3t = mapreduce(c -> sample(m8_3t, sampler, nsamples), chainscat, 1:nchains)

# Results rethinking

m83rethinking = "
      mean   sd  5.5% 94.5% n_eff Rhat
alpha 0.09 1.63 -2.13  2.39   959    1
sigma 2.04 2.05  0.68  4.83  1090    1
";

# End of `08/m8.3t.jl`
