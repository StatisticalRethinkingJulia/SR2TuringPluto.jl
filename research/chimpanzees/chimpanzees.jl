# ### m10.3t.jl

using Pkg, DrWatson

@quickactivate "SR2TuringPluto"
using Turing
using StatisticalRethinking
Turing.setprogress!(false)

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

data = (df.pulled_left, df.condition, df.prosoc_left)
m10_3t = ppl10_3(data...)
nchains = 4; sampler = NUTS(0.65); nsamples=2000
chns10_3t = mapreduce(c -> sample(m10_3t, sampler, nsamples), chainscat, 1:nchains)

# Rethinking result

m_10_03t_result = "
      Mean StdDev lower 0.89 upper 0.89 n_eff Rhat
 a    0.05   0.13      -0.15       0.25  3284    1
 bp   0.62   0.22       0.28       0.98  3032    1
 bpC -0.11   0.26      -0.53       0.29  3184    1
";

chains = sample(m10_3t, NUTS(1000, .65), MCMCThreads(), 1000, 4)

# method for MCMCChains
function pointwise_loglikes(chain::Chains, data, ll_fun)
    samples = Array(Chains(chain, :parameters).value)
    pointwise_loglikes(samples, data, ll_fun)
end

# generic method for arrays
function pointwise_loglikes(samples::Array{Float64,3}, data, ll_fun)
    n_data = length(data)
    n_samples, n_chains = size(samples)[[1,3]]
    pointwise_lls = fill(0.0, n_data, n_samples, n_chains)
    for c in 1:n_chains 
        for s in 1:n_samples
            for d in 1:n_data
                pointwise_lls[d,s,c] = ll_fun(samples[s,:,c], data[d])
            end
        end
    end
    return pointwise_lls
end

function compute_loo(psis_output, pointwise_lls)
    dims = size(pointwise_lls)
    lwp = deepcopy(pointwise_lls)
    lwp += psis_output.weights;
    lwpt = reshape(lwp, dims[1], dims[2] * dims[3])';
    loos = reshape(logsumexp(lwpt; dims=1), size(lwpt, 2));
    return sum(loos)
end

# compute the pointwise log likelihoods where indices correspond to [data, sample, chain]
pointwise_lls = pointwise_loglikes(chains, df.pulled_left, (p,d)->logpdf(Normal(p...), d))

# compute the psis object
psis_output = psis(pointwise_lls)

# return loo based on Rob's example
loo = compute_loo(psis_output, pointwise_lls)

