
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Clip-02-03-06t.jl"

md"## snippet 2.3"

begin
	p_grid = range(0, 1, length = 20)
	prior = ones(20)
	likelihood = pdf.(Binomial.(9, p_grid), 6)
	posterior = likelihood .* prior
	posterior ./= sum(posterior)
end

md"## snippet 2.4"

plot(p_grid, posterior, m = 3, legend = false,
    xlabel = "probability of water", ylabel = "posterior probability", title = "20 points")

md"## snippet 2.5"

begin
	prior2 = ifelse.(p_grid .< 0.5, 0, 1)
	prior3 = @. exp(-5 * abs(p_grid - 0.5))  # @. means broadcast everything that follows
end

md"## snippet 2.6"

@model globethrowing(W, L) = begin
    p ~ Uniform(0, 1)
    W ~ Binomial(W + L, p)
end

m = globethrowing(6, 3)

r = quap(m)

d = MvNormal([r.coef.p], √collect(reshape(r.vcov, 1, 1)))

s = rand(d, 10000)';

histogram(collect(s), normalize = :probability)

plot(p_grid, posterior, m = 3)

md"## End clip-02-03-06t.jl"

