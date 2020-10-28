
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using LinearAlgebra
	using Turing
	using StatisticalRethinking
	Turing.turnprogress(false)
end

md"## Clip-04-26-36t.jl"

md"### snippet 4.26"

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
end;

md"### snippet 4.27 - 4.29"

@model function ppl4_1(heights)
    μ ~ Normal(178, 20)
    σ ~ Uniform(0, 50)
    heights .~ Normal(μ, σ)
end

md"### snippet 4.30"

m4_1t = ppl4_1(df.height)

precis(m4_1t)

md"##### Evaluating above cell might occasionally fail in Optim. Below cell should work."

begin
	start = [mean(Normal(178, 20)), mean(Uniform(0, 50))]
	q4_1t = quap(m4_1t, start)
end

md"### snippet 4.31"

@model function ppl4_2(heights)
    μ ~ Normal(178, 0.1)
    σ ~ Uniform(0, 50)
    heights .~ Normal(μ, σ)
end

begin
	m4_2t = ppl4_2(df.height)
	q4_2t = quap(m4_2t, NelderMead())
end

typeof(q4_2t)

q4_1t.params

md"### snippets 4.32, 4.33"

q4_1t.vcov

diag(q4_1t.vcov)

cov2cor(Matrix(q4_1t.vcov), sqrt.(diag(q4_1t.vcov)))

md"### snippets 4.34 - 4.36"

begin
	res = rand(q4_1t.distr, 10_000)
	quap4_1t = DataFrame(res', ["μ", "σ"])
	Text(precis(quap4_1t; io=String))
end

begin
	chns4_1t = sample(m4_1t, NUTS(), 1000)
end

md"## End of clip-04-26-36t.jl"

