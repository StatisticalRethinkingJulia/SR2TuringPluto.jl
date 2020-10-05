
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

md"## Fig3.3t.jl"

begin
	p_grid = range(0, step=0.000001, stop=1)
	prior = ones(length(p_grid))
	likelihood = pdf.(Binomial.(3, p_grid), 3)
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
	samples = sample(p_grid, Weights(posterior), length(p_grid));
end;

begin
	b1 = quantile(samples, [0.25, 0.75])
	b2 = hpdi(samples, alpha=0.5)

	p1 = plot_density_interval(samples, b1;
	  xlab="Proportion water (p)", title="50% PI");
	p2 = plot_density_interval(samples, b2;
	  xlab="Proportion water (p)", title="50% HPDI");
end;

plot(p1, p2, layout=(1, 2))

md"## End of Fig3.3t.jl"

