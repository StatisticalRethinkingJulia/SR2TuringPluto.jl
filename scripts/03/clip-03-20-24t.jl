
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

md"## Clip-03-20-24t.jl"

md"### snippet 3.20"

pdf.(Binomial(2, 0.7), 0:2)

md"### snippet 3.21"

rand(Binomial(2, 0.7))

md"### snippet 3.22"

rand(Binomial(2, 0.7), 10)

begin
	N = 100000
	d = rand(Binomial(2, 0.7), N);
	[length(filter(e -> e == i, d)) for i in 0:2] / N |> display

	d = rand(Binomial(9, 0.7), N);
	h1 = histogram([filter(e -> e == i, d) for i in 0:9];
	 bins=-0.5:1:9.5, color=:lightblue, leg=false, xticks=0:9, bar_width=0.2)
end

md"## End of clip-03-20-24t.jl"

