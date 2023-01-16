### A Pluto.jl notebook ###
# v0.19.3

html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

using Markdown
using InteractiveUtils

# ╔═╡ 9f79ad1a-9bdf-4fe9-b4e4-94fd01b420d2
using Pkg, DrWatson

# ╔═╡ 7f95b8ef-71ee-44c1-adb4-9bfffcca8246
begin
    using Distributions
    using StatsPlots
    using StatsBase
    using LaTeXStrings
    using CSV
    using DataFrames
    using LinearAlgebra
    using Logging
    using Random
    using Turing
    using Dagitty
    using StatisticalRethinking
	using StatisticalRethinkingPlots
end

# ╔═╡ f10da537-82aa-4826-b7fd-3f8f86e11ec7
md"### Set defaults for plot and logging."

# ╔═╡ 001d0e14-cc62-415a-8a07-90aa69d426c6
begin
	default(label=false)
	Logging.disable_logging(Logging.Warn);
end;

# ╔═╡ d7e8f8b2-63c6-4e2b-a660-da89f55426ac
md"## 5.1 Spurious association."

# ╔═╡ b5557598-d1b1-40e8-892e-23bb8ed8fcbe
md"### Code 5.1"

# ╔═╡ f50a5e4f-05ce-49ad-a3e2-68399c0d690b
begin
	d = CSV.read(sr_datadir("WaffleDivorce.csv"), DataFrame)
	d[!,:D] = standardize(ZScoreTransform, d.Divorce)
	d[!,:M] = standardize(ZScoreTransform, d.Marriage)
	d[!,:A] = standardize(ZScoreTransform, d.MedianAgeMarriage);
end

# ╔═╡ 1bf80d02-910a-4d81-b5be-8a59a43150b3
md"### Code 5.2"

# ╔═╡ fe9ac4a0-4235-4717-ba2f-cd3fa75469d5
std(d.MedianAgeMarriage)

# ╔═╡ c1bf4af0-5785-413c-b9e5-8413bb0d3f76
md"### Code 5.3"

# ╔═╡ f9eec178-10a9-4d1b-984f-3e7803289da4
@model function m5_1(A, D)
    σ ~ Exponential(1)
    a ~ Normal(0, 0.2)
    bA ~ Normal(0, 0.5)
    μ = @. a + bA * A
    D ~ MvNormal(μ, σ)
end

# ╔═╡ 323ed22e-fd8b-499e-bd61-13ad7fd774bb
begin
	m5_1t = sample(m5_1(d.A, d.D), NUTS(), 1000)
	m5_1_df = DataFrame(m5_1t)
	prior = sample(m5_1([0], [0]), Prior(), 1000)
	prior_df = DataFrame(prior)
	PRECIS(prior_df)
end

# ╔═╡ 6b3ea4ab-9129-40e6-ad98-c3c3bb616749
md"### Code 5.4"

# ╔═╡ 41e30acc-3c9a-4524-9189-d45aa23b9f3f
let
	# calculate μ for every prior sample on age=-2 and age=2
	
	bounds = [-2, 2]
	μ = StatisticalRethinking.link(prior_df, [:a, :bA], bounds)
	μ = hcat(μ...);

	p = plot(xlab="Median age marriage (std)", ylab="Divorce rate (std)")
	for μₚ ∈ first(eachrow(μ), 50)
		plot!(bounds, μₚ; c=:black, alpha=0.3)
	end
end

# ╔═╡ a35679e8-d86b-4687-845d-eb2cfd292761
md"### Code 5.5"

# ╔═╡ 98258325-f115-49fc-8129-16633da2a1cb
let
	A_seq = range(-3, 3.2; length=30)

	μ = StatisticalRethinking.link(m5_1_df, [:a, :bA], A_seq)
	μ = hcat(μ...)
	μ_mean = mean.(eachcol(μ))
	μ_PI = PI.(eachcol(μ))
	μ_PI = vcat(μ_PI'...)

	@df d scatter(:A, :D; xlab="Median age marriage (std)",
		ylab="Divorce rate (std)")
	plot!(A_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.2)
end

# ╔═╡ 91f20b09-e6af-4a68-a8b5-b21fa07e5d8e
md"### Code 5.6"

# ╔═╡ beb5142b-e810-4534-bff5-2b6788e7957f
@model function m5_2(M, D)
    σ ~ Exponential(1)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    μ = @. a + bM * M
    D ~ MvNormal(μ, σ)
end

# ╔═╡ e9728aea-1da4-4ee0-a37c-66263625680a
begin
	m5_2t = sample(m5_2(d.M, d.D), NUTS(), 1000)
	m5_2_df = DataFrame(m5_2t)
	PRECIS(m5_2_df)
end

# ╔═╡ d61e3eef-4c40-4eab-bc8a-b4a03f613730
let
	M_seq = range(-1.74, 2.8; length=30)

	μ = StatisticalRethinking.link(m5_2_df, [:a, :bM], M_seq)
	μ = hcat(μ...)
	μ_mean = mean.(eachcol(μ))
	μ_PI = PI.(eachcol(μ))
	μ_PI = vcat(μ_PI'...)

	@df d scatter(:M, :D; xlab="Marriage rate (std)", ylab="Divorce rate (std)")
	plot!(M_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.2)
end

# ╔═╡ 1d1ca39e-5e8a-4902-9391-582e9dda2bd2
md"### Code 5.7"

# ╔═╡ 9b8944d7-4f95-4184-b2a8-eab619932ef3
let
	g = Dagitty.DAG(:A => :M, :A => :D, :M => :D)
	drawdag(g, [0, 1, 2], [0, 1, 0])
end

# ╔═╡ 06e4d0cc-16a5-4703-9b05-7463f3e53a14
md"### Code 5.8"

# ╔═╡ 5accef48-803e-40e0-bffc-572a0ae22d52
let
	g = Dagitty.DAG(:A => :M, :A => :D)
	implied_conditional_independencies(g)
end

# ╔═╡ 673d4f3b-2a7e-4b2d-aa4f-be02fbf10783
md"### Code 5.9"

# ╔═╡ f4d42055-4f3d-4c95-ab43-a5e7c80da6b5
let
	g = Dagitty.DAG(:A => :M, :A => :D, :M => :D)
	implied_conditional_independencies(g)
end

# ╔═╡ 953261e8-8596-4a07-9d34-5de87ad85a00
md"### Code 5.10"

# ╔═╡ 8e8abfd3-5757-4de6-8720-464f04a41178
@model function m5_3(A, M, D)
    σ ~ Exponential(1)
    a ~ Normal(0, 0.2)
    bA ~ Normal(0, 0.5)
    bM ~ Normal(0, 0.5)
    μ = @. a + bA * A + bM * M
    D ~ MvNormal(μ, σ)
end

# ╔═╡ 80f84b76-6027-4b76-955b-0f399097a684
begin
	m5_3t = sample(m5_3(d.A, d.M, d.D), NUTS(), 1000)
	m5_3_df = DataFrame(m5_3t)
	PRECIS(m5_3_df)
end

# ╔═╡ 9bf794d7-f391-4395-a8e5-75d62f6f22b0
md"### Code 5.11"

# ╔═╡ 971d8348-9d68-48aa-b07b-edbcf5dfc209
coeftab_plot(m5_1_df, m5_2_df, m5_3_df; pars=(:bA, :bM),
	names=["m5.1", "m5.2", "m5.3"])

# ╔═╡ fbe4c1a9-6ab4-4bc6-a475-9b8cf5f19fb6
md"### Code 5.12"

# ╔═╡ 869b634b-4047-49a6-a17f-0987f277d6a8
let
	N = 50
	age = rand(Normal(), N)
	mar = rand.(Normal.(-age))
	div = rand.(Normal.(age));

	s1 = DataFrame(sample(m5_1(age, div), NUTS(), 1000))
	s2 = DataFrame(sample(m5_2(mar, div), NUTS(), 1000))
	s3 = DataFrame(sample(m5_3(age, mar, div), NUTS(), 1000));
	coeftab_plot(s1, s2, s3; pars=(:bA, :bM), names=["s1", "s2", "s3"])
end

# ╔═╡ b9313a4f-eb80-4d2e-9d27-59eee5eaa8c7
let
	N = 50
	age = rand(Normal(), N)
	mar = rand.(Normal.(-age))
	div = rand.(Normal.(age .+ mar));

	s1 = DataFrame(sample(m5_1(age, div), NUTS(), 1000))
	s2 = DataFrame(sample(m5_2(mar, div), NUTS(), 1000))
	s3 = DataFrame(sample(m5_3(age, mar, div), NUTS(), 1000));
	coeftab_plot(s1, s2, s3; pars=(:bA, :bM), names=["s1", "s2", "s3"])
end

# ╔═╡ e8460cac-dfd6-4f6c-bba2-1c06bf5e2d42
md"### Code 5.13"

# ╔═╡ ae36190b-665a-4c2c-801f-4fd890ddc8dc
@model function m5_4(A, M)
    σ ~ Exponential(1)
    a ~ Normal(0, 0.2)
    bAM ~ Normal(0, 0.5)
    μ = @. a + bAM * A
    M ~ MvNormal(μ, σ)
end

# ╔═╡ 278a9162-c619-4c5d-aa5c-c785751bdeab
begin
	m5_4t = sample(m5_4(d.A, d.M), NUTS(), 1000)
	m5_4_df = DataFrame(m5_4t);
end

# ╔═╡ 3ca7c4ed-fee3-4beb-b279-6ed00a19e91e
md"### Code 5.14"

# ╔═╡ 2393365d-6a7a-455e-95e1-78dc66796fdd
let
	mu = StatisticalRethinking.link(m5_4_df, [:a, :bAM], d.A);
	mu = hcat(mu...)
	mu_mean = mean.(eachcol(mu))
	mu_resid = mu_mean .- d.M;

	# +
	# Side-note: how to plot the residuals
	# getting yerr - list of 2-tuples with distance to the regression line
	yerr = collect(zip(-clamp.(mu_resid, -Inf, -0.0), clamp.(mu_resid, 0, Inf)));

	plot(d.A, mu_mean; xlab="Age at marriage (std)", ylab="Marriage rate (std)")
	scatter!(d.A, d.M)
	scatter!(d.A, d.M; yerr=yerr, markersize=0)
end

# ╔═╡ 88daa32f-56ea-47d6-8d70-8abca0f2b0d8
md"### Code 5.15"

# ╔═╡ 6500ebc7-55a6-477f-a245-23c30764d275
let
	
	# explicit link form before I improved it
	
	mu = [
		@. r.a + r.bA * d.A + r.bM * d.M
		for r ∈ eachrow(m5_3_df)
	]

	mu = vcat(mu'...)
	mu_mean = mean.(eachcol(mu))
	mu_PI = PI.(eachcol(mu))
	mu_PI = vcat(mu_PI'...);

	D_sim = [
		rand(MvNormal((@. r.a + r.bA * d.A + r.bM * d.M), r.σ))
		for r ∈ eachrow(m5_3_df)
	]
	D_sim = vcat(D_sim'...);
	D_PI = PI.(eachcol(D_sim))
	D_PI = vcat(D_PI'...);
	# -

	# Code 5.16

	yerr = mu_PI[:,2] .- mu_mean
	scatter(d.D, mu_mean; xlab="Observed divorce", ylab="Predicted divorce",
		yerr=yerr)
	plot!(x->x; style=:dash)

	# Code 5.17

	loc_flags = d.Loc .∈ (["ID", "UT", "RI", "ME"],);
	loc_idxes = findall(loc_flags);
	anns = [
		(d.D[idx] - 0.1, mu_mean[idx], (d.Loc[idx], 8))
		for idx in loc_idxes
	]
	annotate!(anns)
end

# ╔═╡ 2184b190-1e08-4aa2-8902-6a4739c64fc1
md"### Code 5.18"

# ╔═╡ 40af7e49-d7a3-4dd9-aeab-c5805f71a94c
let
	N = 100
	x_real = rand(Normal(), N)
	x_spur = rand.(Normal.(x_real))
	y = rand.(Normal.(x_real))
	df = DataFrame(:y => y, :x_real => x_real, :x_spur => x_spur)
	PRECIS(df)
end

# ╔═╡ 6a3f4b38-3640-46e6-b203-0aa25d7be651
md"### Code 5.19"

# ╔═╡ 4f071c8c-92bf-4dc5-b85b-f87ad251a63a
@model function m5_3A(A, M, D)
	# A → D ← M
	σ ~ Exponential(1)
	a ~ Normal(0, 0.2)
	bA ~ Normal(0, 0.5)
	bM ~ Normal(0, 0.5)
	μ = @. a + bA * A + bM * M
	D ~ MvNormal(μ, σ)
	# A → M
	σ_M ~ Exponential(1)
	aM ~ Normal(0, 0.2)
	bAM ~ Normal(0, 0.5)
	μ_M = @. aM + bAM * A
	M ~ MvNormal(μ_M, σ_M)
end

# ╔═╡ 58ef9eed-31b5-4459-9350-ca2a9270b50e
begin
	d1 = CSV.read(sr_datadir("WaffleDivorce.csv"), DataFrame)
	d2 = DataFrame(
		:D => standardize(ZScoreTransform, d1.Divorce),
		:M => standardize(ZScoreTransform, d1.Marriage),
		:A => standardize(ZScoreTransform, d1.MedianAgeMarriage),
	);

	m5_3At = sample(m5_3A(d2.A, d2.M, d2.D), NUTS(), 1000)
	m5_3A_df = DataFrame(m5_3At)
	PRECIS(m5_3A_df)
end

# ╔═╡ 11936398-443c-42c7-ba97-0915afd13ce2
md"### Code 5.20 - 5.22"

# ╔═╡ 80c1a904-0054-47c1-b8cf-1f53e99325ba
let
	A_seq = range(-2, 2; length=30);
	global s_M, s_D = [], []

	for r ∈ eachrow(m5_3A_df)
		M = rand(MvNormal((@. r.aM + r.bAM * A_seq), r.σ_M))
		D = rand(MvNormal((@. r.a + r.bA * A_seq + r.bM * M), r.σ))
		push!(s_M, M)
		push!(s_D, D)
	end

	s_M = vcat(s_M'...)
	s_D = vcat(s_D'...);
	μ_D = mean.(eachcol(s_D))
	PI_D = vcat(PI.(eachcol(s_D))'...)

	plot(
		A_seq, [μ_D, μ_D]; 
		fillrange=PI_D, fillalpha=0.2, color=:black,
		xlab="manupulated A", ylab="counterfactual D",
		title="Total counterfactual effect of A on D"
	)
end

# ╔═╡ 2732c579-c6ad-479c-94c8-16f1dbaf8c6d
md"### Code 5.23"

# ╔═╡ 74d8915b-8da8-4b8d-a2e9-cbba05b86909
let
	sim2_A = @. ([20, 30] - 26.1) / 1.24;
	s2_M, s2_D = [], []

	for r ∈ eachrow(m5_3A_df)
		M = rand(MvNormal((@. r.aM + r.bAM * sim2_A), r.σ_M))
		D = rand(MvNormal((@. r.a + r.bA * sim2_A + r.bM * M), r.σ))
		push!(s2_M, M)
		push!(s2_D, D)
	end

	s2_M = vcat(s2_M'...)
	s2_D = vcat(s2_D'...);
	mean(s2_D[:,2] - s2_D[:,1])
	# -

	# Code 5.24

	# +
	M_seq = range(-2, 2; length=30)
	s_D = []

	for r ∈ eachrow(m5_3A_df)
		# A is zero, so, we drop it from the μ term
		D = rand(MvNormal((@. r.a + r.bM * M_seq), r.σ))
		push!(s_D, D)
	end

	s_D = vcat(s_D'...);

	μ_D = mean.(eachcol(s_D))
	PI_D = vcat(PI.(eachcol(s_D))'...)

	plot(
		M_seq, [μ_D, μ_D]; 
		fillrange=PI_D, fillalpha=0.2, color=:black,
		xlab="manupulated M", ylab="counterfactual D",
		title="Total counterfactual effect of M on D"
	)
end

# ╔═╡ 2a1871be-088f-4fdf-a80b-2d27a752e86c
md"### Code 5.25"

# ╔═╡ f45d0c18-8572-48c3-af57-f9698bb77c3c
A_seq = range(-2, 2; length=30);

# ╔═╡ 66db4fc0-a89b-4011-a781-2591e4256c77
let
	μ_M = mean.(eachcol(s_M))
	PI_M = vcat(PI.(eachcol(s_M))'...)

	plot(
		A_seq, [μ_M, μ_M]; 
		fillrange=PI_M, fillalpha=0.2, color=:black,
		xlab="manupulated A", ylab="counterfactual M",
		title="Total counterfactual effect of A on M"
	)
end

# ╔═╡ 21c81fc5-f7d9-411a-9d5c-01375cc82cb8
md"### Code 5.26"

# ╔═╡ 23bf6317-2baa-49c2-9251-179cbe755158
let
	s_M = [
		rand(MvNormal((@. r.aM + r.bAM * A_seq), r.σ_M))
		for r ∈ eachrow(m5_3A_df)
	]
	s_M = vcat(s_M'...);

	# Code 5.27

	s_D = [
		rand(MvNormal((@. r.a + r.bA * A_seq + r.bM * M), r.σ))
		for (r, M) ∈ zip(eachrow(m5_3A_df), eachrow(s_M))
	]
	s_D = vcat(s_D'...);
end

# ╔═╡ Cell order:
# ╠═9f79ad1a-9bdf-4fe9-b4e4-94fd01b420d2
# ╠═7f95b8ef-71ee-44c1-adb4-9bfffcca8246
# ╟─f10da537-82aa-4826-b7fd-3f8f86e11ec7
# ╠═001d0e14-cc62-415a-8a07-90aa69d426c6
# ╟─d7e8f8b2-63c6-4e2b-a660-da89f55426ac
# ╟─b5557598-d1b1-40e8-892e-23bb8ed8fcbe
# ╠═f50a5e4f-05ce-49ad-a3e2-68399c0d690b
# ╟─1bf80d02-910a-4d81-b5be-8a59a43150b3
# ╠═fe9ac4a0-4235-4717-ba2f-cd3fa75469d5
# ╟─c1bf4af0-5785-413c-b9e5-8413bb0d3f76
# ╠═f9eec178-10a9-4d1b-984f-3e7803289da4
# ╠═323ed22e-fd8b-499e-bd61-13ad7fd774bb
# ╟─6b3ea4ab-9129-40e6-ad98-c3c3bb616749
# ╠═41e30acc-3c9a-4524-9189-d45aa23b9f3f
# ╟─a35679e8-d86b-4687-845d-eb2cfd292761
# ╠═98258325-f115-49fc-8129-16633da2a1cb
# ╟─91f20b09-e6af-4a68-a8b5-b21fa07e5d8e
# ╠═beb5142b-e810-4534-bff5-2b6788e7957f
# ╠═e9728aea-1da4-4ee0-a37c-66263625680a
# ╠═d61e3eef-4c40-4eab-bc8a-b4a03f613730
# ╟─1d1ca39e-5e8a-4902-9391-582e9dda2bd2
# ╠═9b8944d7-4f95-4184-b2a8-eab619932ef3
# ╟─06e4d0cc-16a5-4703-9b05-7463f3e53a14
# ╠═5accef48-803e-40e0-bffc-572a0ae22d52
# ╟─673d4f3b-2a7e-4b2d-aa4f-be02fbf10783
# ╠═f4d42055-4f3d-4c95-ab43-a5e7c80da6b5
# ╟─953261e8-8596-4a07-9d34-5de87ad85a00
# ╠═8e8abfd3-5757-4de6-8720-464f04a41178
# ╠═80f84b76-6027-4b76-955b-0f399097a684
# ╟─9bf794d7-f391-4395-a8e5-75d62f6f22b0
# ╠═971d8348-9d68-48aa-b07b-edbcf5dfc209
# ╟─fbe4c1a9-6ab4-4bc6-a475-9b8cf5f19fb6
# ╠═869b634b-4047-49a6-a17f-0987f277d6a8
# ╠═b9313a4f-eb80-4d2e-9d27-59eee5eaa8c7
# ╟─e8460cac-dfd6-4f6c-bba2-1c06bf5e2d42
# ╠═ae36190b-665a-4c2c-801f-4fd890ddc8dc
# ╠═278a9162-c619-4c5d-aa5c-c785751bdeab
# ╟─3ca7c4ed-fee3-4beb-b279-6ed00a19e91e
# ╠═2393365d-6a7a-455e-95e1-78dc66796fdd
# ╟─88daa32f-56ea-47d6-8d70-8abca0f2b0d8
# ╠═6500ebc7-55a6-477f-a245-23c30764d275
# ╟─2184b190-1e08-4aa2-8902-6a4739c64fc1
# ╠═40af7e49-d7a3-4dd9-aeab-c5805f71a94c
# ╟─6a3f4b38-3640-46e6-b203-0aa25d7be651
# ╠═4f071c8c-92bf-4dc5-b85b-f87ad251a63a
# ╠═58ef9eed-31b5-4459-9350-ca2a9270b50e
# ╟─11936398-443c-42c7-ba97-0915afd13ce2
# ╠═80c1a904-0054-47c1-b8cf-1f53e99325ba
# ╠═66db4fc0-a89b-4011-a781-2591e4256c77
# ╟─2732c579-c6ad-479c-94c8-16f1dbaf8c6d
# ╠═74d8915b-8da8-4b8d-a2e9-cbba05b86909
# ╟─2a1871be-088f-4fdf-a80b-2d27a752e86c
# ╠═f45d0c18-8572-48c3-af57-f9698bb77c3c
# ╟─21c81fc5-f7d9-411a-9d5c-01375cc82cb8
# ╠═23bf6317-2baa-49c2-9251-179cbe755158
