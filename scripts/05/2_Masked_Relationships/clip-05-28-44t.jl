
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Clip-05-28-44.jl"

md"### snippets 5.28 - 5.32"

md"##### Notice that `missing` values are dropped!"

begin
	df = CSV.read(sr_datadir("milk.csv"), DataFrame; missingstring = "NA")
	dropmissing!(df, [:kcal_per_g, :neocortex_perc, :mass])
	df.K = zscore(df.kcal_per_g)
	df.N = zscore(df.neocortex_perc)
	df.M = zscore(log.(df.mass))
end;

md"### snippet 5.33"

@model function m5_5_draft(N, K)
    a ~ Normal(0, 1)
    bN ~ Normal(0, 1)
    σ ~ Exponential(1)
    μ = lin(a, bN, N)
    K ~ MvNormal(μ, σ)
end

begin
	m5_5_draftt = m5_5_draft(df.N, df.K)
	quap5_5_draftt = quap(m5_5_draftt)
end

md"### snippet 5.34"

xseq = -2:0.1:2

begin
	prior = sample(m5_5_draftt, Prior(), 50) |> DataFrame
	μ = lin(prior.a', prior.bN', xseq)

	fig1 = plot(
		xlab="neocortex percent (std)", ylab="kcal per g (std)",
		title="a ~ Normal(0, 1)\nbN ~ Normal(0, 1)",
		legend = false, ylims = extrema(xseq))
	for c in eachcol(μ)
		plot!(xseq, c, color = :black, alpha = 0.3)
	end
	plot!()
end

md"### snippet 5.35"

@model function m5_5(N, K)
    a ~ Normal(0, 0.2)
    bN ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, bN, N)
    K ~ MvNormal(μ, σ)
end

begin
	m5_5t = m5_5(df.N, df.K)
	chns5_5t = sample(m5_5t, NUTS(0.65), 2000)
end

begin
	post = sample(m5_5_draftt, NUTS(0.65), 200) |> DataFrame
	μ_post = lin(post.a', post.bN', xseq)

	fig2 = plot(
		xlab="neocortex percent (std)", ylab="kcal per g (std)",
		title="a ~ Normal(0, 0.2)\nbN ~ Normal(0, 0.5)",
		legend = false, ylims = extrema(xseq))
	for c in eachcol(μ_post)
		plot!(xseq, c, color = :black, alpha = 0.3)
	end
	plot(fig1, fig2, layout=(1, 2))
end

md"### snippet 5.36"

begin
	post5_5t = DataFrame(Array(chns5_5t), names(chns5_5t, [:parameters]))
	Text(precis(post5_5t; io=String))
end

begin
	quap5_5t = quap(m5_5t)
	Text(precis(quap5_5t; io=String))
end

md"### snippet 5.37"

begin
	post1 = DataFrame(rand(quap5_5t.distr, 1000)', quap5_5t.params)
	μ1 = lin(post1.a', post1.bN', xseq) |> meanlowerupper

	fig3 = plot(xlab="neocortex percent (std)", ylab="kcal per gram (std)")
	scatter!(df.N, df.K, legend = false)
	plot!(xseq, μ1.mean, ribbon = (μ1.mean .- μ1.lower, μ1.upper .- μ1.mean))
end

md"### snippet 5.38"

@model function m5_6(M, K)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, bM, M)
    K ~ MvNormal(μ, σ)
end

begin
	m5_6t = m5_6(df.M, df.K)
	Text(precis(quap(m5_6t); io=String))
end

begin
	quap5_6t = quap(m5_6t)
	post2 = DataFrame(rand(quap5_6t.distr, 1000)', quap5_6t.params)
	μ2 = lin(post2.a', post2.bM', xseq) |> meanlowerupper

	fig4 = plot(xlab="log body mass (std)", ylab="kcal per gram (std)")
	scatter!(df.N, df.K, legend = false)
	plot!(xseq, μ2.mean, ribbon = (μ2.mean .- μ2.lower, μ2.upper .- μ2.mean))
end


md"### snippet 5.39"

@model function m5_7(M, N, K)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    bN ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, bM, M, bN, N)
    K ~ MvNormal(μ, σ)
end

begin
	m5_7t = m5_7(df.M, df.N, df.K)
	quap5_7t = quap(m5_7t)
	post5_7t = DataFrame(rand(quap5_7t.distr, 1000)', quap5_7t.params)
	Text(precis(post5_7t; io=String))
end

md"### snippet 5.40"

begin
	quap_res, fig = plotcoef([m5_5t, m5_6t, m5_7t], [:a, :bN, :bM, :σ])
	plot(fig)
end

quap_res

md"### snippet 5.41"

begin
	mu2 = lin(post5_7t.a', post5_7t.bN', xseq, post5_7t.bM', zeros(length(xseq))) |> meanlowerupper
	fig5 = plot(xlab="neocortex percent (std)", ylab="kcal per gram (std)",
		title="Counterfactual (M=0)", leg=false)
	plot!(xseq, mu2.mean, ribbon = (mu2.upper .- mu2.mean, mu2.mean .- mu2.lower))

	mu3 = lin(post5_7t.a', post5_7t.bM', xseq, post5_7t.bN', zeros(length(xseq))) |> meanlowerupper
	fig6 = plot(xlab="neocortex percent (std)", ylab="log body mass (std)", 
		title="Counterfactual (N=0)", leg=false)
	plot!(xseq, mu3.mean, ribbon = (mu3.upper .- mu3.mean, mu3.mean .- mu3.lower))
	
	plot(fig3, fig4, fig5, fig6, layout=(2,2))
end

md"### snippet 5.42"

n = 100

begin
	# M -> K <- N
	# M -> N
	M1 = randn(n)
	N1 = rand.(Normal.(M1))
	K1 = rand.(Normal.(N1 .- M1))
	d_sim1 = DataFrame((; K1, N1, M1))
	Text(precis(d_sim1; io=String))
end

md"### snippet 5.43"

begin
	# M -> K <- N
	# N -> M
	N2 = randn(n)
	M2 = rand.(Normal.(N2))
	K2 = rand.(Normal.(N2 .- M2))
	d_sim2 = DataFrame((; K2, N2, M2))
	Text(precis(d_sim2; io=String))
end

begin
	# M -> K <- N
	# M -> U <- N
	U3 = randn(n)
	N3 = rand.(Normal.(U3))
	M3 = rand.(Normal.(U3))
	K3 = rand.(Normal.(N3 .- M3))
	d_sim3 = DataFrame((; K3, N3, M3))
	Text(precis(d_sim3; io=String))
end

md"### snippet 5.44"

md"##### TODO working with DAGs."

md"## End of clip-05-28-44.jl"

