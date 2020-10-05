
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end;

md"## Clip-04-48-49t.jl"

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)     # either
end;

@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* (weights .- x̄)
    heights ~ MvNormal(μ, σ)
end

md"### snippets 4.48, 4.49"

begin
	figs = Vector{Plots.Plot{Plots.GRBackend}}(undef, 3)
		for (indx, N) in enumerate([50, 150, 352])

		dN = df[1:N, :]
		quapN = quap(m4_3(dN.weight, dN.height))

		post = rand(quapN.distr, 20)
		post = DataFrame(post', quapN.params)

		figs[indx] = scatter(dN.weight, dN.height)
		for p in eachrow(post)
			plot!(x, p.a .+ p.b .* (x .- mean(dN.weight)), color = "black", alpha = 0.3)
		end
		plot!(legend = false, xlabel = "weight", ylabel = "height")
	end
	plot(figs..., layout=(1, 3))
end

md"## End of clip-04-48-49t.jl"

