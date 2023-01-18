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

# ╔═╡ 80deafdd-ce8d-4b7d-9651-7ba281b38873
using Pkg, DrWatson

# ╔═╡ 12d8d9bd-c72e-43e3-a183-6611619e774b
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

# ╔═╡ 39cb3b35-d8c1-4dd7-a4f1-83fc2fe20a64
md"#### Set defaults for plotting and logging."

# ╔═╡ 9b65e1ff-e7eb-4800-9862-718c167a8e32
begin
	default(label=false)
	Logging.disable_logging(Logging.Warn);
end

# ╔═╡ 6b9df7a9-556b-4b73-b67f-2ac8c48aa0c7
md"## 5.2 Masked relationship."

# ╔═╡ 2e0252f4-58c2-4fa1-a063-c2b5dd6a274a
md"### Code 5.28"

# ╔═╡ 3de73995-159f-4a62-8ca3-91a7c5c378e0
begin
	d = CSV.read(sr_datadir("milk.csv"),  missingstring="NA", DataFrame)

	# get rid of dots in column names

		rename!(n -> replace(n, "." => "_"), d)

	describe(d)

end

# ╔═╡ 84025c5a-49d6-49bb-ada7-1a33e50fe0c4
md"### Code 5.29"

# ╔═╡ 7674fdcb-c7e8-4b87-ba55-ea05cd1a1f91
begin
	d[!,:K] = standardize(ZScoreTransform, d.kcal_per_g)
	d[!,:M] = standardize(ZScoreTransform, log.(d.mass))

	# column contains missing values, need to propagate them on standartization
	
	d[!,:N] = d.neocortex_perc
	non_miss = findall(!ismissing, d.N);
	d[non_miss,:N] = standardize(ZScoreTransform, disallowmissing(d.N[non_miss]));
end;

# ╔═╡ adac6445-d51d-43a4-89fd-19a14bde89c4
md"### Code 5.30"

# ╔═╡ 37f91d86-0f12-4b87-935c-83e0bdb61dd0
@model function model_m5_5_draft(N, K)
    a ~ Normal(0, 1)
    bN ~ Normal(0, 1)
    σ ~ Exponential(1)
    μ = @. a + bN * N
    K ~ MvNormal(μ, σ)
end

# ╔═╡ 01a4d9e6-493e-47bf-af7b-e1e826475d52
try
    m5_5_draft = sample(model_m5_5_draft(d.N, d.K), NUTS(), 1000)
catch e
    if isa(e, MethodError)
        s = sprint(showerror, e)
        println(s)
    end
end

# ╔═╡ d1be5297-6a50-4028-8927-698ee6c3e804
md"### Code 5.31"

# ╔═╡ 9099b078-1ae1-4b78-bbc7-c9d5e7d9fea9
d.neocortex_perc

# ╔═╡ dbe69ae2-0646-4370-9d3e-918fc2db686c
md"### Code 5.32"

# ╔═╡ db826c7f-7396-4350-ac24-2aeddf50ba66
dcc = d[completecases(d[!,[:K,:N,:M]]),:];

# ╔═╡ d6134b92-8951-49dd-a448-6092f2ecabbe
md"### Code 5.33"

# ╔═╡ 8d10f16e-748c-4cfd-9297-f1d2eff99d13
m5_5_draft = sample(model_m5_5_draft(dcc.N, dcc.K), NUTS(), 1000);

# ╔═╡ e540f095-d956-4ff9-9f5b-f6ec07d5cad5
md"### Code 5.34"

# ╔═╡ 399e259d-6997-4017-98c2-a4e09df1acfb
let
	prior = sample(model_m5_5_draft(dcc.N, dcc.K), Prior(), 1000)
	prior_df = DataFrame(prior)
	xseq = [-2, 2]
	μ = StatisticalRethinking.link(prior_df, [:a, :bN], xseq)
	μ = hcat(μ...);

	p = plot(; xlim=xseq, ylim=xseq, 
		xlab="neocortex percent (std)", ylab="kilocal per g (std)", 
		title=L"a \sim \mathcal{N}(0,1), bN \sim \mathcal{N}(0,1)"
	)
	for y ∈ first(eachrow(μ), 50)
		plot!(p, xseq, y; c=:black, alpha=0.3)
	end
	p
end

# ╔═╡ 1b9660fb-4074-459c-ab07-ffcad4d45a26
md"### Code 5.35"

# ╔═╡ 4793fe6a-735b-4949-8e27-a1a8f6c8e3c8
@model function model_m5_5(N, K)
    a ~ Normal(0, 0.2)
    bN ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = @. a + bN * N
    K ~ MvNormal(μ, σ)
end

# ╔═╡ ea7a7a99-29d5-4589-81bd-88d8aa64d885
begin
	m5_5 = sample(model_m5_5(dcc.N, dcc.K), NUTS(), 1000)
	m5_5_df = DataFrame(m5_5)
	PRECIS(m5_5_df)
end

# ╔═╡ 34b193e3-a4c3-40a9-b913-ef980e0aee35
let
	prior = sample(model_m5_5(dcc.N, dcc.K), Prior(), 1000)
	prior_df = DataFrame(prior)
	xseq = [-2, 2]


	μ = StatisticalRethinking.link(prior_df, [:a, :bN], xseq)
	μ = hcat(μ...);

	p2 = plot(; xlim=xseq, ylim=xseq, 
		xlab="neocortex percent (std)", ylab="kilocal per g (std)", 
		title=L"a \sim \mathcal{N}(0,0.2), bN \sim \mathcal{N}(0,0.5)"
	)
	for y ∈ first(eachrow(μ), 50)
		plot!(p2, xseq, y; c=:black, alpha=0.3)
	end
end

# ╔═╡ 2b52e28e-57b1-4f94-9672-9ad1ca528594
md"### Code 5.36"

# ╔═╡ e292be02-0e44-4a6a-8de0-4046eac7c8b3
PRECIS(m5_5_df)

# ╔═╡ f48c5a95-b9e1-4783-943d-fbd51f8ccc79
md"### Code 5.37"

# ╔═╡ 075ce6cf-2b2f-4d88-97a4-d08a6d28bb02
let
	xseq = range(minimum(dcc.N) - 0.15, maximum(dcc.N) + 0.15; length=30)
	μ = StatisticalRethinking.link(m5_5_df, [:a, :bN], xseq);
	μ = hcat(μ...)
	μ_mean = mean.(eachcol(μ))
	μ_PI = PI.(eachcol(μ))
	μ_PI = vcat(μ_PI'...)

	@df dcc scatter(:N, :K; xlab="neocortex percent (std)",
		ylab="kilocal per g (std)")
	plot!(xseq, [μ_mean, μ_mean]; lw=2, fillrange=μ_PI, fillalpha=0.2, color=:black)
end

# ╔═╡ 74cbc4ab-b2d8-4e8e-8496-89b2c6e47ac7
md"### Code 5.38"

# ╔═╡ e51a1fa2-8716-49e1-806d-c54610212711
@model function model_m5_6(M, K)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = @. a + bM * M
    K ~ MvNormal(μ, σ)
end

# ╔═╡ c49709e2-4893-422f-8936-615590d7773a
begin
	m5_6 = sample(model_m5_6(dcc.M, dcc.K), NUTS(), 1000)
	m5_6_df = DataFrame(m5_6)
	PRECIS(m5_6_df)
end

# ╔═╡ 17eea9d5-802c-49fd-b29f-d311811044a1
let
	xseq = range(minimum(dcc.M) - 0.15, maximum(dcc.M) + 0.15; length=30)
	μ = StatisticalRethinking.link(m5_6_df, [:a, :bM], xseq);
	μ = hcat(μ...)
	μ_mean = mean.(eachcol(μ))
	μ_PI = PI.(eachcol(μ))
	μ_PI = vcat(μ_PI'...)

	@df dcc scatter(:M, :K; xlab="log body mass (std)", ylab="kilocal per g (std)")
	plot!(xseq, [μ_mean, μ_mean]; lw=2, fillrange=μ_PI, fillalpha=0.2, color=:black)
end

# ╔═╡ a1babcb0-d1a3-454f-95a4-7d7091e98f3e
@model function model_m5_7(N, M, K)
    a ~ Normal(0, 0.2)
    bN ~ Normal(0, 0.5)
    bM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = @. a + bN * N + bM * M
    K ~ MvNormal(μ, σ)
end

# ╔═╡ fed67c7b-23f5-4e95-9bd7-cf174946b1e1
begin
	m5_7 = sample(model_m5_7(dcc.N, dcc.M, dcc.K), NUTS(), 1000)
	m5_7_df = DataFrame(m5_7)
	PRECIS(m5_7_df)
end

# ╔═╡ f4dd2f6b-0991-45cc-82ec-d4c361e8067b
md"### Code 5.40"

# ╔═╡ bc6c18e4-ee46-4f0a-9ee3-c1ea7351184e
coeftab_plot(m5_7_df, m5_6_df, m5_5_df; pars=(:bM, :bN),
	names=("m5.7", "m5.6", "m5.5"))

# ╔═╡ 845129b4-db30-4d58-b9ad-682f5a9b15e3
md"### Code 5.41"

# ╔═╡ d8e83c21-fbeb-40d3-aa80-b810271d7399
md"##### The code in the book corresponds to the bottom-right figure, which keeps N=0 (despite stated in the text)."

# ╔═╡ 51ec6d48-9608-4021-9a3a-b89a4e3555d5
md"##### Below is the code to produce the bottom-left figure (M=0)."

# ╔═╡ a2f99116-de8d-4183-8725-9a22081d0b79
let
	xseq = range(minimum(dcc.N) - 0.15, maximum(dcc.N) + 0.15; length=30)
	μ = StatisticalRethinking.link(m5_7_df, [:a, :bN], xseq);
	μ = hcat(μ...)
	μ_mean = mean.(eachcol(μ))
	μ_PI = PI.(eachcol(μ))
	μ_PI = vcat(μ_PI'...)

	plot(title="Counterfactual holding M=0", 
		xlab="neocortex percent (std)", ylab="kilocal per g (std)")
	plot!(xseq, [μ_mean, μ_mean]; lw=2, fillrange=μ_PI, fillalpha=0.2, color=:black)

	# +
	xseq = range(minimum(dcc.M) - 0.15, maximum(dcc.M) + 0.15; length=30)
	μ = StatisticalRethinking.link(m5_7_df, [:a, :bM], xseq);
	μ = hcat(μ...)
	μ_mean = mean.(eachcol(μ))
	μ_PI = PI.(eachcol(μ))
	μ_PI = vcat(μ_PI'...)

	plot(title="Counterfactual holding N=0", 
		xlab="log body mass (std)", ylab="kilocal per g (std)")
	plot!(xseq, [μ_mean, μ_mean]; lw=2, fillrange=μ_PI, fillalpha=0.2, color=:black)
end

# ╔═╡ 343e3408-595e-47ba-8100-a99cb4dc0781
md"### Code 5.42"

# ╔═╡ 4d6dc4aa-fb17-44fb-8029-4775866be39d
let
	# M → K ← N
	# M → N
	
	n = 100
	M = rand(Normal(), n)
	N = [rand(Normal(μ)) for μ ∈ M]
	K = [rand(Normal(μ)) for μ ∈ N .- M] 
	d_sim = DataFrame(:K => K, :N => N, :M => M);

	s5 = sample(model_m5_5(d_sim.N, d_sim.K), NUTS(), 1000)
	s6 = sample(model_m5_6(d_sim.M, d_sim.K), NUTS(), 1000)
	s7 = sample(model_m5_7(d_sim.N, d_sim.M, d_sim.K), NUTS(), 1000)
	s5_df = DataFrame(s5)
	s6_df = DataFrame(s6)
	s7_df = DataFrame(s7)
	coeftab_plot(s7_df, s6_df, s5_df; pars=(:bM, :bN), names=("s7", "s6", "s5"))
end

# ╔═╡ 5e9a7797-a8cb-4f30-8982-08f04d251797
md"### Code 5.43"

# ╔═╡ df429c0b-7297-4a80-bfc0-1e68c7112961
let
	# M → K ← N
	# N → M
	n = 100
	N = rand(Normal(), n)
	M = [rand(Normal(μ)) for μ ∈ N]
	K = [rand(Normal(μ)) for μ ∈ N .- M] 
	d_sim2 = DataFrame(:K => K, :N => N, :M => M);

	# M → K ← N
	# M ← U → N
	n = 100
	U = rand(Normal(), n)
	N = [rand(Normal(μ)) for μ ∈ U]
	M = [rand(Normal(μ)) for μ ∈ U]
	K = [rand(Normal(μ)) for μ ∈ N .- M] 
	global d_sim3 = DataFrame(:K => K, :N => N, :M => M);
	# -

	s5 = sample(model_m5_5(d_sim2.N, d_sim2.K), NUTS(), 1000)
	s6 = sample(model_m5_6(d_sim2.M, d_sim2.K), NUTS(), 1000)
	s7 = sample(model_m5_7(d_sim2.N, d_sim2.M, d_sim2.K), NUTS(), 1000)
	s5_df = DataFrame(s5)
	s6_df = DataFrame(s6)
	s7_df = DataFrame(s7)
	coeftab_plot(s7_df, s6_df, s5_df; pars=(:bM, :bN), names=("s7", "s6", "s5"))
end

# ╔═╡ 48ad9b7b-1f32-4acb-9e9c-9f823e1d6255
let
s5 = sample(model_m5_5(d_sim3.N, d_sim3.K), NUTS(), 1000)
	s6 = sample(model_m5_6(d_sim3.M, d_sim3.K), NUTS(), 1000)
	s7 = sample(model_m5_7(d_sim3.N, d_sim3.M, d_sim3.K), NUTS(), 1000)
	s5_df = DataFrame(s5)
	s6_df = DataFrame(s6)
	s7_df = DataFrame(s7)
	coeftab_plot(s7_df, s6_df, s5_df; pars=(:bM, :bN), names=("s7", "s6", "s5"))
end

# ╔═╡ 13d4f85d-1a4e-4209-850d-2ee2ca0890a7
md"### Code 5.44"

# ╔═╡ f3bb3a8d-265e-4053-9917-922805c9d8ad
begin
	dag5_7 = Dagitty.DAG(:M => :K, :N => :K, :M => :N)
	drawdag(dag5_7, [1, 0, 1], [0, 0, 1])
end

# ╔═╡ 83d1ddf2-e9a2-4c13-864e-f9c3978b161f
md"##### EquivalentDAGs is TODO in Dagitty.jl."

# ╔═╡ Cell order:
# ╠═80deafdd-ce8d-4b7d-9651-7ba281b38873
# ╠═12d8d9bd-c72e-43e3-a183-6611619e774b
# ╟─39cb3b35-d8c1-4dd7-a4f1-83fc2fe20a64
# ╠═9b65e1ff-e7eb-4800-9862-718c167a8e32
# ╠═6b9df7a9-556b-4b73-b67f-2ac8c48aa0c7
# ╠═2e0252f4-58c2-4fa1-a063-c2b5dd6a274a
# ╠═3de73995-159f-4a62-8ca3-91a7c5c378e0
# ╟─84025c5a-49d6-49bb-ada7-1a33e50fe0c4
# ╠═7674fdcb-c7e8-4b87-ba55-ea05cd1a1f91
# ╟─adac6445-d51d-43a4-89fd-19a14bde89c4
# ╠═37f91d86-0f12-4b87-935c-83e0bdb61dd0
# ╠═01a4d9e6-493e-47bf-af7b-e1e826475d52
# ╟─d1be5297-6a50-4028-8927-698ee6c3e804
# ╠═9099b078-1ae1-4b78-bbc7-c9d5e7d9fea9
# ╟─dbe69ae2-0646-4370-9d3e-918fc2db686c
# ╠═db826c7f-7396-4350-ac24-2aeddf50ba66
# ╟─d6134b92-8951-49dd-a448-6092f2ecabbe
# ╠═8d10f16e-748c-4cfd-9297-f1d2eff99d13
# ╟─e540f095-d956-4ff9-9f5b-f6ec07d5cad5
# ╠═399e259d-6997-4017-98c2-a4e09df1acfb
# ╟─1b9660fb-4074-459c-ab07-ffcad4d45a26
# ╠═4793fe6a-735b-4949-8e27-a1a8f6c8e3c8
# ╠═ea7a7a99-29d5-4589-81bd-88d8aa64d885
# ╠═34b193e3-a4c3-40a9-b913-ef980e0aee35
# ╠═2b52e28e-57b1-4f94-9672-9ad1ca528594
# ╠═e292be02-0e44-4a6a-8de0-4046eac7c8b3
# ╠═f48c5a95-b9e1-4783-943d-fbd51f8ccc79
# ╠═075ce6cf-2b2f-4d88-97a4-d08a6d28bb02
# ╟─74cbc4ab-b2d8-4e8e-8496-89b2c6e47ac7
# ╠═e51a1fa2-8716-49e1-806d-c54610212711
# ╠═c49709e2-4893-422f-8936-615590d7773a
# ╠═17eea9d5-802c-49fd-b29f-d311811044a1
# ╠═a1babcb0-d1a3-454f-95a4-7d7091e98f3e
# ╠═fed67c7b-23f5-4e95-9bd7-cf174946b1e1
# ╟─f4dd2f6b-0991-45cc-82ec-d4c361e8067b
# ╠═bc6c18e4-ee46-4f0a-9ee3-c1ea7351184e
# ╟─845129b4-db30-4d58-b9ad-682f5a9b15e3
# ╟─d8e83c21-fbeb-40d3-aa80-b810271d7399
# ╟─51ec6d48-9608-4021-9a3a-b89a4e3555d5
# ╠═a2f99116-de8d-4183-8725-9a22081d0b79
# ╟─343e3408-595e-47ba-8100-a99cb4dc0781
# ╠═4d6dc4aa-fb17-44fb-8029-4775866be39d
# ╟─5e9a7797-a8cb-4f30-8982-08f04d251797
# ╠═df429c0b-7297-4a80-bfc0-1e68c7112961
# ╠═48ad9b7b-1f32-4acb-9e9c-9f823e1d6255
# ╠═13d4f85d-1a4e-4209-850d-2ee2ca0890a7
# ╠═f3bb3a8d-265e-4053-9917-922805c9d8ad
# ╟─83d1ddf2-e9a2-4c13-864e-f9c3978b161f
