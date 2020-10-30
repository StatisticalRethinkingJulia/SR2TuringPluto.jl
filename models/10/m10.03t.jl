# ### m10.3t.jl

using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking
Turing.turnprogress(false)

delim=';'
df = CSV.read(sr_datadir("chimpanzees.csv"), DataFrame; delim);

# pulled_left, condition, prosoc_left
@model function ppl10_3(y, x₁, x₂)
    α ~ Normal(0, 10)
    βp ~ Normal(0, 10)
    βpC ~ Normal(0, 10)

    logits = α .+ (βp .+ βpC * x₁) .* x₂
    y .~ BinomialLogit.(1, logits)
end;

m10_3t = ppl10_3(df.pulled_left, df.condition, df.prosoc_left)
nchains = 4; sampler = NUTS(0.65); nsamples=2000
chns10_3t = mapreduce(c -> sample(m10_3t, sampler, nsamples), chainscat, 1:nchains)

# Rethinking result

m_10_03t_result = "
      Mean StdDev lower 0.89 upper 0.89 n_eff Rhat
 a    0.05   0.13      -0.15       0.25  3284    1
 bp   0.62   0.22       0.28       0.98  3032    1
 bpC -0.11   0.26      -0.53       0.29  3184    1
";

# End of m10.03t.jl
