
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

md"## Fig3.2t.jl"

begin
	p_grid = range(0, stop=1, length=10000)
	prior = ones(length(p_grid))
	likelihood = pdf.(Binomial.(9, p_grid), 6)
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
	samples = sample(p_grid, Weights(posterior), length(p_grid));
end;

begin
	N = 10000
	samples2 = sample(p_grid, Weights(posterior), N)
end;

begin
	b1 = mapreduce(p -> p < 0.5 ? 1 : 0, +, samples2) / N
	b2 = mapreduce(p -> (p > 0.5 && p < 0.75) ? 1 : 0, +, samples2) / N
	b3 = quantile(samples2, 0.8)
	b4 = quantile(samples2, [0.1, 0.9])

	p1 = plot_density_interval(samples2, [0.0, 0.5],
	  xlab="Proportion water (p)");
	p2 = plot_density_interval(samples2, [0.5, 0.75],
	  xlab="Proportion water (p)");
	p3 = plot_density_interval(samples2, [0.0, b3], 
	  xlab="Proportion water (p)");
	p4 = plot_density_interval(samples2, b4, 
	  xlab="Proportion water (p)")
end;

plot(p1, p2, p3, p4, layout=(2, 2))

md"## End of Fig3.2t.jl"

