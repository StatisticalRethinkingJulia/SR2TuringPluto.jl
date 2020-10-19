
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Clip-05-50-54t.jl"

md"### snippet 5.50"

begin
	df = CSV.read(sr_datadir("milk.csv"), DataFrame; missingstring = "NA")
	df.clade |> unique |> sort
end;

md"### snippet 5.51"

md"##### You can't just turn Strings into Integers in Julia but hashing them should give the same result"

df.clade_id = Int.(indexin(df.clade, unique(df.clade)))

md"### snippet 5.52"

df.K = zscore(df.kcal_per_g)

md"##### Clade DynamicPPL model"

@model function ppl5_9(clade_id, K)
    σ ~ Exponential(1)
    α ~ filldist(Normal(0, 0.5), length(unique(clade_id)))
    μ = α[clade_id]
    K ~ MvNormal(μ, σ)
end

begin
	m5_9t = ppl5_9(df.clade_id, df.K)
	quap5_9t = quap(m5_9t)
	post5_9t = DataFrame(rand(quap5_9t.distr, 1000)', quap5_9t.params)
	Text(precis(post5_9t; io=String))
end


md"### snippet 5.54"

df.house = rand(1:4, nrow(df));

@model function ppl5_10(clade_id, house, K)
    σ ~ Exponential(1)
    a ~ filldist(Normal(0, 0.5), length(unique(clade_id)))
    h ~ filldist(Normal(0, 0.5), length(unique(house)))
    μ = a[clade_id] + h[house]
    K ~ MvNormal(μ, σ)
end

begin
	m5_10t = ppl5_10(df.clade_id, df.house, df.K)
	quap5_10t = quap(m5_10t)
	Text(precis(DataFrame(rand(quap5_10t.distr, 1000)', quap5_10t.params); io=String))
end

md"## End of clip-05-45-54.jl"

