
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Fig 2.5.2t"

md"##### This clip is only intended to illustrate a part of Fig 2.5."

md"##### Create a Turing modelt:"

@model globe_toss(W, L) = begin
    p ~ Uniform(0, 1)
    W ~ Binomial(W + L, p)
end

m = globe_toss(6, 3);

r = quap(m)

d = Normal(r.coef.p, √collect(reshape(r.vcov, 1, 1))[1])

s = rand(d, 10000);

h1 = histogram(s, normalize = :pdf, lab="sample density")

begin
	k = 3
	n = 6
	chn = sample(globe_toss(n, k), NUTS(0.65), 11000)
end

plot(chn; seriestype=:traceplot)

begin
	histogram(s, normalize = :pdf, lab="sample density")
	plot!(chn; seriestype=:density, lab="Chn density")
end

md"## End of Fig2.5.2t.jl"

