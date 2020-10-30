# ### m10.xxt.jl

using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking
Turing.turnprogress(false)

# outcome and predictor almost perfectly associated

x = repeat([-1], 9); append!(x, repeat([1],11))
y = repeat([0], 10); append!(y, repeat([1],10))

@model ppl10_xx(x,y) = begin
    α ~ Normal(0,10)
    β ~ Normal(0,10)

    logits = α .+ β * x

    y .~ BinomialLogit.(1, logits)
end

m10_xxt = ppl10_xx(x,y)
nchains = 4; sampler = NUTS(0.65); nsamples=2000
chns10_xxt = mapreduce(c -> sample(m10_xxt, sampler, nsamples), chainscat, 1:nchains)

# Stan results

m_10_x,_results = "
    mean   sd   5.5% 94.5% n_eff Rhat
 a -5.09 4.08 -12.62 -0.25   100 1.00
 b  7.86 4.09   2.96 15.75   104 1.01
";

# End of 10/m10.xxt.jl
