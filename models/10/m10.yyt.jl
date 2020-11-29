# ### m10.yyt.jl

using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking
Turing.setprogress!(false)

delim=';'
df = CSV.read(sr_datadir("UCBadmit.csv"), DataFrame; delim);

@model ppl10_yy(admit, reject) = begin
   α₁ ~ Normal(0,100)
   α₂ ~ Normal(0,100)

   for i ∈ 1:length(admit)
       λₐ = exp(α₁)
       λᵣ = exp(α₂)
       admit[i] ~ Poisson(λₐ)
       reject[i] ~ Poisson(λᵣ)
   end
end;

m10_yyt = ppl10_yy(df.admit, df.reject)
nchains = 4; sampler = NUTS(0.65); nsamples=2000
chns10_yyt = mapreduce(c -> sample(m10_yyt, sampler, nsamples), chainscat, 1:nchains)

# Rethinking/CmdStan result

m_10_yyt_result = "
    mean   sd 5.5% 94.5% n_eff Rhat
 a1 4.99 0.02 4.95  5.02  2201    1
 a2 5.44 0.02 5.41  5.47  2468    1
";

# End of 10/m10.yyt.jl
