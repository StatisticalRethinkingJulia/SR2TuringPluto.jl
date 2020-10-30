# ### m8.1t

using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking
Turing.turnprogress(false)

@model ppl8_2(y) = begin

    σ ~ FlatPos(0.0) # improper prior with probability one everywhere above 0.0
    α ~ Flat() # improper prior with pobability one everywhere

    y .~ Normal(α, σ)
end

y = [-1,1];

# Sample

m8_2t = ppl8_2(y)
nchains = 4; sampler = NUTS(0.65); nsamples=2000
chns8_2t = mapreduce(c -> sample(m8_2t, sampler, nsamples), chainscat, 1:nchains)

# End of `08/m8.2t.jl`
