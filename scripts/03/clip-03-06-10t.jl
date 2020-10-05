
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

md"## Clip-03-06-10t.jl"

md"### snippet 3.2"

begin
	p_grid = range(0, step=0.001, stop=1)
	prior = ones(length(p_grid))
	likelihood = pdf.(Binomial.(9, p_grid), 6)
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
end;

md"### snippet 3.3"

md"##### Draw 10000 samples from this posterior distribution."

begin
	N = 10000
	samples = sample(p_grid, Weights(posterior), N)
end;

md"##### Store samples in an MCMCChains.Chains object."

chn = MCMCChains.Chains(reshape(samples, N, 1, 1), [:p]);

md"##### Describe the chain."

chn

md"##### Plot the chain."

plot(chn)

md"### snippet 3.6"

v1 = sum(posterior[filter(i -> p_grid[i] < 0.5, 1:length(p_grid))])

md"### snippet 3.7"

b1 = mapreduce(p -> p < 0.5 ? 1 : 0, +, samples) / N

md"### snippet 3.8"

b2 = mapreduce(p -> (p > 0.5 && p < 0.75) ? 1 : 0, +, samples) / N

md"### snippet 3.9"

b3 = quantile(samples, 0.8)

md"### snippet 3.10"

b4 = quantile(samples, [0.1, 0.9])

begin
	p1 = plot_density_interval(samples, [0.0, 0.5],
	  xlab="Proportion water (p)");
	p2 = plot_density_interval(samples, [0.5, 0.75],
	  xlab="Proportion water (p)");
	p3 = plot_density_interval(samples, [0.0, b3], 
	  xlab="Proportion water (p)");
	p4 = plot_density_interval(samples, b4, 
	  xlab="Proportion water (p)");
	plot(p1, p2, p3, p4, layout=(2, 2))
end

md"## End of clip-03-06-10t.jl"

