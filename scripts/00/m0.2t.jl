
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
	using StatsPlots
end

using MLJ, MLJModels

md"## m0.2t.jl"

md"##### This notebook takes a look at the [rstar()](https://arxiv.org/pdf/2003.07900.pdf) diagnostic. Currently I can't get MLJ & MLJModels to work (serialization issue?)."

md"##### Define a simple Normal model with unknown mean and variance."

@model gdemo(x, y) = begin
  s ~ InverseGamma(2, 3)
  m ~ Normal(0, sqrt(s))
  x ~ Normal(m, sqrt(s))
  y ~ Normal(m, sqrt(s))
end

md"#####  Run sampler, collect results."

chns = mapreduce(c -> sample(gdemo(1.5, 2), NUTS(0.65), 2000), chainscat, 1:4)

md"##### Plot the chains."

plot(chns; seriestype=:traceplot)

plot(chns; seriestype=:density)

md"##### Below will likely fail!"

chn ... # sampling results of multiple chains

md"##### Select classifier used to compute the diagnostic."

classif = @load XGBoostClassifier

md"##### stimate diagnostic."

Rs = rstar(classif, chn)

R = mean(Rs)

md"##### Visualize distribution."

histogram(Rs)

md"## End m0.2t.jl"

