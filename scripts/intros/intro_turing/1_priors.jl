
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing 
	using StatisticalRethinking
	Turing.turnprogress(false);
end;

md"## Intro clip 1_priors.jl"

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)
	Text(precis(df; io=String))
end

@model function m4_3(weights, heights)
    a ~ Normal(mean(weights), 50)
    b ~ Normal(0, 5)
    σ ~ Uniform(0, 20)
    μ = a .+ b .* weights
    for i in eachindex(heights)
        heights[i] ~ Normal(μ[i], σ)
    end
end

begin
	m4_3t = m4_3(df.weight, df.height)
	priors4_3t = sample(m4_3t, Prior(), 50)
end

begin
	priors4_3t_df = DataFrame(Array(priors4_3t, :parameters), names(priors4_3t, :parameters))
	Text(precis(priors4_3t_df; io=String))
end

begin
	plot(xlims=(25, 70), ylims=(-10, 300), title="Possible prior regression lines for this model", leg=false)
	for row in eachrow(priors4_3t_df)
		plot!(df.weight, row.a .+ row.b * df.weight)
	end
	scatter!(df.weight, df.height)
end

md"## End of intro clip 1_priors.jl"

