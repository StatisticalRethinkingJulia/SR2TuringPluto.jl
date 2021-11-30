### A Pluto.jl notebook ###
# v0.12.3

using Markdown
using InteractiveUtils

# ╔═╡ 92be606a-0998-11eb-053e-cb42246e0598
using Pkg, DrWatson

# ╔═╡ 92be94e0-0998-11eb-03fc-4b2293b5e32a
begin
	@quickactivate "SR2TuringPluto"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ 8cc034b8-0998-11eb-2c71-578b14536d5d
md"## Clip-05-28-44.jl"

# ╔═╡ 92c90f38-0998-11eb-086f-4b4900891b2f
md"### snippets 5.28 - 5.32"

# ╔═╡ 94042790-099b-11eb-1a63-6dc1af913def
md"##### Notice that `missing` values are dropped!"

# ╔═╡ 92d87388-0998-11eb-1881-11abc4118c34
begin
	df = CSV.read(sr_datadir("milk.csv"), DataFrame; missingstring = "NA")
	dropmissing!(df, [:kcal_per_g, :neocortex_perc, :mass])
	df.K = zscore(df.kcal_per_g)
	df.N = zscore(df.neocortex_perc)
	df.M = zscore(log.(df.mass))
end;

# ╔═╡ 55142898-0999-11eb-1cff-6d82e287b64a
md"### snippet 5.33"

# ╔═╡ fdde27d4-099a-11eb-2156-9b15bba277f2
@model function m5_5_draft(N, K)
    a ~ Normal(0, 1)
    bN ~ Normal(0, 1)
    σ ~ Exponential(1)
    μ = lin(a, bN, N)
    K ~ MvNormal(μ, σ)
end

# ╔═╡ fdde6780-099a-11eb-0fd8-b597b7b34d15
begin
	m5_5_draftt = m5_5_draft(df.N, df.K)
	quap5_5_draftt = quap(m5_5_draftt)
end

# ╔═╡ e17b5820-099d-11eb-242b-c1f8c9e2fc14
md"### snippet 5.34"

# ╔═╡ d03a085c-09a2-11eb-3ba3-d18810ca8ab5
xseq = -2:0.1:2

# ╔═╡ e4d38704-099d-11eb-0c30-9d7f0556f5b1
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

# ╔═╡ e4d3b814-099d-11eb-2d2d-bbf106ef33b6
md"### snippet 5.35"

# ╔═╡ e4d468b6-099d-11eb-2b83-3b2d265a4fad
@model function m5_5(N, K)
    a ~ Normal(0, 0.2)
    bN ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, bN, N)
    K ~ MvNormal(μ, σ)
end

# ╔═╡ e4ea98f4-099d-11eb-0ca5-bd416c74dbab
begin
	m5_5t = m5_5(df.N, df.K)
	chns5_5t = sample(m5_5t, NUTS(0.65), 2000)
end

# ╔═╡ 82cf498e-09a3-11eb-1edf-4f32d126c20e
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

# ╔═╡ e4ebf03e-099d-11eb-16e3-93fea8ed8f8a
md"### snippet 5.36"

# ╔═╡ e4f60d74-099d-11eb-324e-31bef9d5181f
begin
	post5_5t = DataFrame(Array(chns5_5t), names(chns5_5t, [:parameters]))
	Text(precis(post5_5t; io=String))
end

# ╔═╡ 26000e2c-099f-11eb-05cc-bb0bede678d1
begin
	quap5_5t = quap(m5_5t)
	Text(precis(quap5_5t; io=String))
end

# ╔═╡ 4010b8b4-09a1-11eb-173a-2bb20afc882d
md"### snippet 5.37"

# ╔═╡ 44a91bfa-09a1-11eb-3bd9-af4110eab915
begin
	post1 = DataFrame(rand(quap5_5t.distr, 1000)', quap5_5t.params)
	μ1 = lin(post1.a', post1.bN', xseq) |> meanlowerupper

	fig3 = plot(xlab="neocortex percent (std)", ylab="kcal per gram (std)")
	scatter!(df.N, df.K, legend = false)
	plot!(xseq, μ1.mean, ribbon = (μ1.mean .- μ1.lower, μ1.upper .- μ1.mean))
end

# ╔═╡ 44a9500c-09a1-11eb-36e3-750c47541694
md"### snippet 5.38"

# ╔═╡ 44aa2266-09a1-11eb-3d17-5d3f2da14256
@model function m5_6(M, K)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, bM, M)
    K ~ MvNormal(μ, σ)
end

# ╔═╡ 44bf965a-09a1-11eb-2d3b-95b5e58e1edd
begin
	m5_6t = m5_6(df.M, df.K)
	Text(precis(quap(m5_6t); io=String))
end

# ╔═╡ 753a9688-09a4-11eb-1bce-654aeb39c454
begin
	quap5_6t = quap(m5_6t)
	post2 = DataFrame(rand(quap5_6t.distr, 1000)', quap5_6t.params)
	μ2 = lin(post2.a', post2.bM', xseq) |> meanlowerupper

	fig4 = plot(xlab="log body mass (std)", ylab="kcal per gram (std)")
	scatter!(df.N, df.K, legend = false)
	plot!(xseq, μ2.mean, ribbon = (μ2.mean .- μ2.lower, μ2.upper .- μ2.mean))
end


# ╔═╡ 44c06a94-09a1-11eb-243e-e94c4ea7b657
md"### snippet 5.39"

# ╔═╡ 44d48ee8-09a1-11eb-0dd6-c5aa78fc84ec
@model function m5_7(M, N, K)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    bN ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, bM, M, bN, N)
    K ~ MvNormal(μ, σ)
end

# ╔═╡ 44d79f14-09a1-11eb-2122-916596c22454
begin
	m5_7t = m5_7(df.M, df.N, df.K)
	quap5_7t = quap(m5_7t)
	post5_7t = DataFrame(rand(quap5_7t.distr, 1000)', quap5_7t.params)
	Text(precis(post5_7t; io=String))
end

# ╔═╡ 2aa59598-09a2-11eb-0854-5b7b1a9aa1b2
md"### snippet 5.40"

# ╔═╡ 99bfbc84-09a1-11eb-01a6-35c26e589d7b
begin
	quap_res, fig = plotcoef([m5_5t, m5_6t, m5_7t], [:a, :bN, :bM, :σ])
	plot(fig)
end

# ╔═╡ 1982c47a-09a2-11eb-357c-1b16ba7aacb1
quap_res

# ╔═╡ 59c2eef4-09a5-11eb-2a39-b75a6de102ad
md"### snippet 5.41"

# ╔═╡ 5e7f6ba2-09a5-11eb-1168-7bbd0f34e93b
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

# ╔═╡ 405a080e-0995-11eb-01c0-3ff35fe8c925
md"### snippet 5.42"

# ╔═╡ 99aa57b4-09a7-11eb-2080-01a888bd0999
n = 100

# ╔═╡ 99aa96fa-09a7-11eb-17c6-bb516e055371
begin
	# M -> K <- N
	# M -> N
	M1 = randn(n)
	N1 = rand.(Normal.(M1))
	K1 = rand.(Normal.(N1 .- M1))
	d_sim1 = DataFrame((; K1, N1, M1))
	Text(precis(d_sim1; io=String))
end

# ╔═╡ 99b7614e-09a7-11eb-21b4-716994e5df9d
md"### snippet 5.43"

# ╔═╡ 99ba1a76-09a7-11eb-1729-339be9f12454
begin
	# M -> K <- N
	# N -> M
	N2 = randn(n)
	M2 = rand.(Normal.(N2))
	K2 = rand.(Normal.(N2 .- M2))
	d_sim2 = DataFrame((; K2, N2, M2))
	Text(precis(d_sim2; io=String))
end

# ╔═╡ 99d1457a-09a7-11eb-18c8-bb41b89014f3
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

# ╔═╡ 99d510ba-09a7-11eb-0fcf-97d4242bb925
md"### snippet 5.44"

# ╔═╡ 99e1a37a-09a7-11eb-3839-87157105995e
md"##### TODO working with DAGs."

# ╔═╡ 99eecdfe-09a7-11eb-12c8-fd7d822403c1
md"## End of clip-05-28-44.jl"

# ╔═╡ Cell order:
# ╟─8cc034b8-0998-11eb-2c71-578b14536d5d
# ╠═92be606a-0998-11eb-053e-cb42246e0598
# ╠═92be94e0-0998-11eb-03fc-4b2293b5e32a
# ╟─92c90f38-0998-11eb-086f-4b4900891b2f
# ╟─94042790-099b-11eb-1a63-6dc1af913def
# ╠═92d87388-0998-11eb-1881-11abc4118c34
# ╟─55142898-0999-11eb-1cff-6d82e287b64a
# ╠═fdde27d4-099a-11eb-2156-9b15bba277f2
# ╠═fdde6780-099a-11eb-0fd8-b597b7b34d15
# ╟─e17b5820-099d-11eb-242b-c1f8c9e2fc14
# ╠═d03a085c-09a2-11eb-3ba3-d18810ca8ab5
# ╠═e4d38704-099d-11eb-0c30-9d7f0556f5b1
# ╟─e4d3b814-099d-11eb-2d2d-bbf106ef33b6
# ╠═e4d468b6-099d-11eb-2b83-3b2d265a4fad
# ╠═e4ea98f4-099d-11eb-0ca5-bd416c74dbab
# ╠═82cf498e-09a3-11eb-1edf-4f32d126c20e
# ╟─e4ebf03e-099d-11eb-16e3-93fea8ed8f8a
# ╠═e4f60d74-099d-11eb-324e-31bef9d5181f
# ╠═26000e2c-099f-11eb-05cc-bb0bede678d1
# ╟─4010b8b4-09a1-11eb-173a-2bb20afc882d
# ╠═44a91bfa-09a1-11eb-3bd9-af4110eab915
# ╟─44a9500c-09a1-11eb-36e3-750c47541694
# ╠═44aa2266-09a1-11eb-3d17-5d3f2da14256
# ╠═44bf965a-09a1-11eb-2d3b-95b5e58e1edd
# ╠═753a9688-09a4-11eb-1bce-654aeb39c454
# ╟─44c06a94-09a1-11eb-243e-e94c4ea7b657
# ╠═44d48ee8-09a1-11eb-0dd6-c5aa78fc84ec
# ╠═44d79f14-09a1-11eb-2122-916596c22454
# ╟─2aa59598-09a2-11eb-0854-5b7b1a9aa1b2
# ╠═99bfbc84-09a1-11eb-01a6-35c26e589d7b
# ╠═1982c47a-09a2-11eb-357c-1b16ba7aacb1
# ╟─59c2eef4-09a5-11eb-2a39-b75a6de102ad
# ╠═5e7f6ba2-09a5-11eb-1168-7bbd0f34e93b
# ╟─405a080e-0995-11eb-01c0-3ff35fe8c925
# ╠═99aa57b4-09a7-11eb-2080-01a888bd0999
# ╠═99aa96fa-09a7-11eb-17c6-bb516e055371
# ╟─99b7614e-09a7-11eb-21b4-716994e5df9d
# ╠═99ba1a76-09a7-11eb-1729-339be9f12454
# ╠═99d1457a-09a7-11eb-18c8-bb41b89014f3
# ╟─99d510ba-09a7-11eb-0fcf-97d4242bb925
# ╟─99e1a37a-09a7-11eb-3839-87157105995e
# ╟─99eecdfe-09a7-11eb-12c8-fd7d822403c1
