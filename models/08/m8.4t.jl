# ### m8.4t

using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking
Turing.turnprogress(false)

# Can't really set a Uniform[-Inf,Inf] on σ 

# Turing model
@model ppl8_4(y) = begin
    α₁ ~ Uniform(-3000, 1000)
    α₂ ~ Uniform(-1000, 3000)
    σ ~ truncated(Cauchy(0,1), 0, Inf)

    y .~ Normal(α₁ + α₂, σ)
end

# Observations

y = rand(Normal(0,1), 100);

# Sample

m8_4t = ppl8_4(y)
nchains = 4; sampler = NUTS(0.65); nsamples=2000
chns8_4t = mapreduce(c -> sample(m8_4t, sampler, nsamples), chainscat, 1:nchains)

# Results rethinking

m8_4rethinking = "
         mean      sd     5.5%   94.5% n_eff Rhat
 a1    -861.15 558.17 -1841.89  -31.04     7 1.43
 a2     861.26 558.17    31.31 1842.00     7 1.43
 sigma    0.97   0.07     0.89    1.09     9 1.17
";

# Describe the posterior samples

chns |> display

# End of `08/m8.4t.jl`
