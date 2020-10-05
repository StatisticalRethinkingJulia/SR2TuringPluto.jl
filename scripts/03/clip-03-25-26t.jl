
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

md"## Clip-03-25-26t.jl"

md"### snippet 3.25"

begin
	N = 10000
	w = rand(Binomial(9, 0.6), N);
	h1 = histogram(w; normalize=:probability, leg=false)
end

md"##### Generate the samples."

begin
	p_grid = range(0, step=0.001, stop=1)
	prior = ones(length(p_grid))
	likelihood = [pdf(Binomial(9, p), 6) for p in p_grid]
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
	samples = sample(p_grid, Weights(posterior), N)
end;

md"### snippet 3.26"

begin
	d = rand.(Binomial.(9, samples));
	h2 = histogram(d; normalize=:probability, 
	  bins=-0.5:1:9.5, leg=false, xticks=0:9, bar_width=0.2)
end

plot(h1, h2, layout=(1, 2))

md"## End of clip-03-25-26t.jl"

