
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using LinearAlgebra
	using Turing
	using StatisticalRethinking
end

md"## Clip-04-26-36t.jl"

md"### snippet 4.26"

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
end;

md"### snippet 4.27 - 4.29"

@model function m4_1(heights)
    μ ~ Normal(178, 20)
    σ ~ Uniform(0, 50)
    heights .~ Normal(μ, σ)
end

md"### snippet 4.30"

begin
	m4_1t = m4_1(df.height)
	quap4_1t = quap(m4_1t)
	# precis(m4_1t)
end

md"##### Evaluating above cell might occasionally fail in Optim. Below cell should work."

begin
	start = [mean(Normal(178, 20)), mean(Uniform(0, 50))]
	quap(m4_1t, start)
	# precis(m4_1t)
end

md"### snippet 4.31"

@model function height(heights)
    μ ~ Normal(178, 0.1)
    σ ~ Uniform(0, 50)
    heights .~ Normal(μ, σ)
end

begin
	m4_2t = height(df.height)
	quap4_2t = quap(m4_2t, NelderMead())
	# precis(m4_2)
end

md"### snippets 4.32, 4.33"

quap4_1t.vcov

diag(quap4_1t.vcov)

cov2cor(Matrix(quap4_1t.vcov), sqrt.(diag(quap4_1t.vcov)))

md"### snippets 4.34 - 4.36"

begin
	post = rand(quap4_1t.distr, 10_000)
	post = DataFrame(post', ["μ", "σ"])
	Text(precis(post; io=String))
end

md"## End of clip-04-26-36t.jl"

