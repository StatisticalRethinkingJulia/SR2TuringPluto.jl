# m2.1t.jl

using Pkg, DrWatson
@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking

Turing.turnprogress(false)

# Define the data

k = 6; n = 9;

# Define the model

@model function ppl2_1(n, k)
  θ ~ Beta(1, 1) # prior
  k ~ Binomial(n, θ) # model
  return k, θ
end;

# Use Turing mcmc

m2_1t = ppl2_1(n, k)
nchains = 4; sampler = NUTS(0.65); nsamples=2000
chns2_1t = mapreduce(c -> sample(m2_1t, sampler, nsamples), chainscat, 1:nchains)
part2_1t = Particles(chns2_1t[:θ])

# End of m2.1t.jl
