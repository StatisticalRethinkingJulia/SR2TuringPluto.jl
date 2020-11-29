# ### m10.10t.jl

using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking
Turing.setprogress!(false)

delim=';'
df = CSV.read(sr_datadir("Kline.csv"), DataFrame; delim);

# New col log_pop, set log() for population data
df[!, :log_pop] = map(x -> log(x), df[:, :population]);

# New col contact_high, set binary values 1/0 if high/low contact
df[!, :contact_high] = map(x -> ifelse(x=="high", 1, 0), df[:, :contact]);

# This is supposed to be a "bad" model since we take non-centered data for the
# predictor log_pop
@model ppl10_10(total_tools, log_pop, contact_high) = begin
    α ~ Normal(0, 100)
    βp ~ Normal(0, 1)
    βc ~ Normal(0, 1)
    βpc ~ Normal(0, 1)

    for i ∈ 1:length(total_tools)
        λ = exp(α + βp*log_pop[i] + βc*contact_high[i] +
            βpc*contact_high[i]*log_pop[i])
        total_tools[i] ~ Poisson(λ)
    end
end;

m10_10t = ppl10_10(df.total_tools, df.log_pop, df.contact_high)
nchains = 4; sampler = NUTS(0.65); nsamples=2000
chns10_10t = mapreduce(c -> sample(m10_10t, sampler, nsamples), chainscat, 1:nchains)

# Rethinking result

m_10_10s_result = "
     mean   sd  5.5% 94.5% n_eff Rhat
 a    0.94 0.37  0.36  1.53  3379    1
 bp   0.26 0.04  0.21  0.32  3347    1
 bc  -0.08 0.84 -1.41  1.23  2950    1
 bpc  0.04 0.09 -0.10  0.19  2907    1
";

# End of m10.10t.jl
