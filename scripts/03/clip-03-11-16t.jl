
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

md"## Clip-03-11-16t.jl"

md"### snippet 3.11"

begin
	p_grid = range(0, step=0.001, stop=1)
	prior = ones(length(p_grid))
	likelihood = [pdf(Binomial(3, p), 3) for p in p_grid]
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
end;

md"##### Draw 10000 samples from this posterior distribution."

begin
	N = 10000
	samples = sample(p_grid, Weights(posterior), N);
end;

md"### snippet 3.13"

hpdi(samples, alpha=0.11)

md"### snippet 3.14"

mode(samples)

md"### snippet 3.15"

mean(samples)

md"### snippet 3.16"

median(samples)

md"##### Plot density."

begin
	density(samples, lab="density")
	vline!(hpdi(samples, alpha=0.5), line=:dash, lab="hpdi")
	vline!(quantile(samples, [0.25, 0.75]), line=:dash, lab="quantile (pi)")
end

md"## End of clip-03-11-16t.jl"

