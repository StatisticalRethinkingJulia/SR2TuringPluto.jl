
using Markdown
using InteractiveUtils

using Pkg, DrWatson

using StatisticalRethinking

md"## Clip-03-02-05t.jl"

@quickactivate "StatisticalRethinkingTuring"

md"### snippet 3.2"

begin
	p_grid = range(0, stop=1, length=1000)
	prior = ones(length(p_grid))
	likelihood = pdf.(Binomial.(9, p_grid), 6)
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
end

md"### snippet 3.3"

md"##### Draw 10000 samples from this posterior distribution."

begin
	N = 10000
	samples = sample(p_grid, Weights(posterior), N)
end;

md"##### Create an MCMCChains.Chains object."

chn = MCMCChains.Chains(reshape(samples, N, 1, 1), [:p]);

md"##### Describe the chain."

chn

md"##### Plot the chain."

p1 = plot(chn; seriestype=:traceplot)

md"### snippet 3.4"

p2 = plot(chn; seriestype=:density)

md"### snippet 3.5"

md"##### Compare with analytical (conjugate) solution."

begin
	w = 6
	n = 9
	x = 0:0.01:1
	plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), lab="Conjugate solution")
	density!(samples, lab="Sample density")
end

begin
	density(samples, lab="Sample2 density")
	vline!(hpdi(samples), lab="hpdi samples2")
	vline!(quantile(samples, [0.25, 0.75]), lab="quantiles [0.25, 0.75]")
end

md"## End of clip-03-02-05t.jl"

