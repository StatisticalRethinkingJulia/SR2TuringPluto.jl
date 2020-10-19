
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Clip-05-45-49.jl"

md"### snippet 5.45"

df = CSV.read(sr_datadir("Howell1.csv"), DataFrame);

md"### snippet 5.46"

begin
	μ_female = rand(Normal(178, 20), 10_000)
	μ_male = rand(Normal(178, 20), 10_000) .+ rand(Normal(0, 10), 10_000)
	Text(precis(DataFrame(; μ_female, μ_male); io=String))
end

md"### snippet 5.47"

df.sex = ifelse.(df.male .== 1, 2, 1);

md"### snippet 5.48"

@model function m5_8(sex, height)
    σ ~ Uniform(0, 50)
    a ~ filldist(Normal(178, 20), 2)
    μ = a[sex]
    height ~ MvNormal(μ, σ)
end

begin
	m5_8t = m5_8(df.sex, df.height)
	quap5_8t = quap(m5_8t, NelderMead())
	Text(precis(quap5_8t))
end

md"### snippet 5.49"

begin
	post5_8t = DataFrame(rand(quap5_8t.distr, 1000)', quap5_8t.params)
	post5_8t.diff_fm = post5_8t[:, "a[1]"] .- post5_8t[:, "a[2]"]
	Text(precis(post5_8t; io=String))
end

md"## End of clip-05-44-49t.jl"

