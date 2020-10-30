# ### m11.5t.jl

using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking
Turing.turnprogress(false)

delim=';'
df = CSV.read(sr_datadir("UCBadmit.csv"), DataFrame; delim);

# Turing model

@model ppl11_5(admit, applications) = begin
    N=length(applications)
    θ ~ truncated(Exponential(1), 0, Inf)
    α ~ Normal(0,2)

    for i ∈ 1:N
        prob = logistic(α)

        # alpha and beta for the BetaBinomial must be provided.
        # The two parameterizations are related by
        # alpha = prob * theta, and beta = (1-prob) * theta.
        # See https://github.com/rmcelreath/rethinking/blob/master/man/dbetabinom.Rd
        alpha = prob * θ
        beta = (1 - prob) * θ

        admit[i] ~ BetaBinomial(applications[i], alpha, beta)
    end
end

# Sample

m11_5t = ppl11_5(df.admit, df.applications)
nchains = 4; sampler = NUTS(0.65); nsamples=2000
chns11_5t = mapreduce(c -> sample(m11_5t, sampler, nsamples), chainscat, 1:nchains)

# Result rethinking

m11_5s = "
         mean   sd  5.5% 94.5% n_eff Rhat
θ        2.74 0.96  1.43  4.37  3583    1
a       -0.37 0.31 -0.87  0.12  3210    1
";

# End of `11/m811.5t.jl`
