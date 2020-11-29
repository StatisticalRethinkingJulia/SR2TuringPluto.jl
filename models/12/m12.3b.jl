# Example model ZIPoisson provided by Chris Fisher

using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking

import Distributions: logpdf, rand

struct ZIPoisson{T1,T2} <: DiscreteUnivariateDistribution
    logλ::T1
    w::T2
end

function logpdf(d::ZIPoisson, y::Int)
    LL = 0.0
    if y == 0
        LLs = zeros(typeof(d.logλ), 2)
        LLs[1] = log(d.w)
        LLs[2] = log(1 - d.w) - exp(d.logλ)
        LL = logsumexp(LLs)
    else
        LL = log(1 - d.w) + logpdf(LogPoisson(d.logλ), y)
    end
    return LL
end

function rand(d::ZIPoisson)
    return rand() <= d.w ? 0 : rand(Poisson(exp(d.logλ)))
end

rand(d::ZIPoisson, N::Int) = map(_->rand(d), 1:N)


inv_logit(x) = 1/(1 + exp(-x))

@model ppl12_3bt(data) = begin
    a1 ~ Normal(1, .5)
    logλ ~ Normal(-1.5, 1)
    w = inv_logit(a1)
    data .~ ZIPoisson(logλ, w)
end

#Random.seed!(74591)
data = rand(ZIPoisson(-1.5, .70), 1000)

n_samples = 3000
n_adapt = 1500
config = NUTS(n_adapt, .65)
n_chains = 4
chns12_3bt = sample(ppl12_3bt(data), config, MCMCThreads(), n_samples, n_chains, progress=true)
chns12_3bt |> display
