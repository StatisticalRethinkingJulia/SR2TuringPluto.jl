
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Clip-04-37-41t.jl"


md"### snippet 4.26"

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)
end;

md"### snippet 4.37"

scatter(df.weight, df.height; lab="Observations", leg=:bottomright)

md"### snippet 4.38"

begin
	N = 100
	a1 = rand(Normal(178.0, 20.0), N)
	b1 = rand(Normal(0.0, 10.0), N)
end;

md"### snippet 4.39"

begin
	plot(xlims = extrema(df.weight), ylims = (-100, 400), xlable = "weight", ylabel = "height")
	hline!([0.0, 272.0])
	for (a, b) in zip(a1, b1)
		plot!(x, a .+ b .* (x .- x̄), color = "black", alpha = 0.2)
	end
	plot!(legend = false)
end

md"### snippet 4.40"

begin
	b = rand(LogNormal(0.0, 1.0), 10_000)
	density(b, xlims = (0, 5))
end

md"### snippet 4.41"

begin
	a2 = rand(Normal(178.0, 20.0), N)
	b2 = rand(LogNormal(0.0, 1.0), N)
	plot(xlims = extrema(df.weight), ylims = (-100, 400), xlable = "weight", ylabel = "height")
	hline!([0.0, 272.0])
	foreach(zip(a2, b2)) do (a, b)
		plot!(x, a .+ b .* (x .- x̄), color = "black", alpha = 0.2)
	end
	plot!(legend = false)
end

md"## End of clip-04-37-41t.jl"

