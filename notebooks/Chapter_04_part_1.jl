### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 2d3d462e-d6bf-4ebd-ab34-db38f4261a11
using Pkg, DrWatson

# ╔═╡ f09ccff8-4036-4322-90e6-cdceb35eadfe
Pkg.activate(expanduser("~/.julia/dev/SR2TuringPluto"))

# ╔═╡ 18b5c138-9324-4c5a-86ca-be1e825ed11e
begin
	using Distributions
	#using Optim
	using StatsPlots
	using StatsBase
	using LaTeXStrings
	using CSV
	using DataFrames
	using StatisticalRethinking
	using LinearAlgebra
	using Logging
	using Random
	using Turing
end

# ╔═╡ 76a7417d-fc80-44c9-b23d-fb75b0e9abba
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 5%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ a2d950f8-d780-4abe-af37-b189bd4a5e88
md"#### Setting default attributes for plots."

# ╔═╡ bbb4a889-77e2-449b-9258-883e7434e641
begin
	default(label=false)
	Logging.disable_logging(Logging.Warn);
end;

# ╔═╡ 0fb4f0c1-df0a-4734-99c9-9a4671d75bb0
md"## 4.1 Why normal distributions are normal."

# ╔═╡ 5fc53e45-ff63-4542-909a-afd5e2716346
md"### Code 4.1"

# ╔═╡ e3709463-a936-485b-8173-0f9fda249417
let
	n = rand(Uniform(-1, 1), 1000, 16);
	pos = sum.(eachrow(n));
	density(pos)
end

# ╔═╡ 58a41999-ae14-423f-b6dd-21fc766374db
md"### Code 4.2"

# ╔═╡ 2f69e9b7-8775-44f9-a493-88d02f8283c7
prod(1 .+ rand(Uniform(0, .1), 12))

# ╔═╡ d694159e-7384-465f-b1ab-23ac6b018323
md"### Code 4.3"

# ╔═╡ 08ee5548-ab90-4eba-b151-a37ce6849df6
let
	u = Uniform(0, .1)
	growth = prod.(eachrow(1 .+ rand(u, 10000, 12)));

	density(growth; label="density")
	
	# Overlay normal distribution
	
	μ = mean(growth)
	σ = std(growth)
	plot!(Normal(μ, σ); label="Normal")
end

# ╔═╡ b5f41800-d83a-4f68-80fc-3b52b9b887ea
md"### Code 4.4"

# ╔═╡ 7fa6c286-638d-4a7d-8fef-1dfa9974b67a
begin
	big = prod.(eachrow(1 .+ rand(Uniform(0, 0.5), 10000, 12)));
	small = prod.(eachrow(1 .+ rand(Uniform(0, 0.01), 10000, 12)));
	density([big, small]; layout=(2, 1))
end

# ╔═╡ 44b3cfbb-cbb9-4640-84f6-097f52a5b117
md"### Code 4.5"

# ╔═╡ b70c459c-8bfc-4ece-9e46-e00ca880be89
density(log.(big))

# ╔═╡ 6f8b69ef-acb5-4db4-b2ae-c08f0fe59e6b
md"## 4.2 A language for describing models."

# ╔═╡ e9862e8d-047b-4d25-992b-2ed5c9dc26bc
md"### Code 4.6"

# ╔═╡ 4d2231d0-9fa9-4600-b325-ff31da8e327f
let
	w = 6
	n = 9
	p_grid = range(0, 1; length=100)
	bin_dens = [pdf(Binomial(n, p), w) for p in p_grid]
	uni_dens = [pdf(Uniform(0, 1), p) for p in p_grid];
	posterior = bin_dens .* uni_dens
	posterior /= sum(posterior);
end

# ╔═╡ d34ed551-5cc3-468c-8c17-35f5c92b2280
md"## 4.3 Gaussian model of height."

# ╔═╡ 687abeee-b831-4e48-aa8a-9a330bc8d5d1
md"### Code 4.7"

# ╔═╡ 496b257c-ea62-40b2-b1c9-e92fdda1d688
begin
	howell1 = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	first(howell1, 5)
end

# ╔═╡ 4753ae08-1b15-4345-8ae6-a4f876457f93
md"### Code 4.8"

# ╔═╡ 50773db0-910c-42ec-b77a-74ad32f23b55
describe(howell1)

# ╔═╡ d117f398-6eba-418a-9696-9f34b820f4b5
md"### Code 4.9"

# ╔═╡ 3dd509b8-4e49-42fc-b3bc-5b2aa8c27542
describe(howell1)

# ╔═╡ 48c8de02-bd59-4bca-b103-d4b770b322d7
md"### Code 4.10"

# ╔═╡ ec60e901-8ea8-4779-bb5d-03d50138e9e6
howell1.height

# ╔═╡ 5737278f-8834-4769-9880-1a8ca580f87c
md"### Code 4.11"

# ╔═╡ 90292dea-4a40-47e6-8269-65f76dd1e50a
begin
	howell1_1 = howell1[howell1.age .>= 18,:]
	describe(howell1_1)
end

# ╔═╡ 7ce59f8a-3f0d-4562-85dd-911a75d68922
md"### Code 4.12"

# ╔═╡ fc7fede8-be68-410a-a147-26ce7ab7e008
plot(Normal(178, 20); xlim=(100, 250))

# ╔═╡ 82b5c59d-0277-4212-af32-1f05459e2623
md"### Code 4.13"

# ╔═╡ ec8fb7f4-83e0-4f3f-8ffa-ec024cfadb88
plot(Uniform(0, 50), xlim=(-10, 60), ylim=(0, 0.03))

# ╔═╡ 83ae8a58-78af-4044-b476-48bff83f818d
md"### Code 4.14"

# ╔═╡ a7100c3c-6fab-43b1-af33-a17780f309e6
let
	size = 10_000
	sample_μ = rand(Normal(178, 20), size)
	sample_σ = rand(Uniform(0, 50), size);
	prior_h = [rand(Normal(μ, σ)) for (μ, σ) in zip(sample_μ, sample_σ)];

	p1 = density(sample_μ; title="μ")
	p2 = density(sample_σ; title="σ")
	p3 = density(prior_h; title="prior_h")

	plot(p1, p2, p3, layout=@layout [a b; c])
end

# ╔═╡ da3e0f8e-9eab-4ebd-bdd7-3738f380fe72
md"### Code 4.15"

# ╔═╡ dfb5d4fe-d810-4705-b03c-c013b8098315
let
	size = 10_000
	sample_μ = rand(Normal(178, 20), size)
	sample_σ = rand(Uniform(0, 50), size);

	sample_μ = rand(Normal(178, 100), size)
	prior_h = [rand(Normal(μ, σ)) for (μ, σ) in zip(sample_μ, sample_σ)];

	density(prior_h)
	vline!([0, 272])
end

# ╔═╡ b72c1949-c95f-4a51-91cb-d47c7cedd4cc
md"### Code 4.16"

# ╔═╡ 98ad03fb-0d07-4cf3-9681-2943a8b4e1fd
begin
	μ_list = range(150, 160; length=100)
	σ_list = range(7, 9; length=100)

	log_likelihood = [
		sum(logpdf.(Normal(μ, σ), howell1_1.height))
		for μ ∈ μ_list, σ ∈ σ_list
	]
	log_prod = log_likelihood .+ [
		logpdf.(Normal(178, 20), μ) + logpdf.(Uniform(0, 50), σ)
		for μ ∈ μ_list, σ ∈ σ_list
	];

	max_prod = maximum(log_prod)
	prob = @. exp(log_prod - max_prod)
	contour(μ_list, σ_list, prob')
end

# ╔═╡ 4c904223-defa-4b45-b8b0-e6cf6f9b7cf7
md"### Code 4.17"

# ╔═╡ 8fcc7434-a202-4e42-9473-a66987c3dd18
md"#### Note the transposition, that's due to Julia matrix order."

# ╔═╡ 56d515ea-74a6-439a-ab34-be0409a3524e
md"### Code 4.18"

# ╔═╡ 12d0c754-ddc9-4ce5-afba-2b07e57fb826
heatmap(μ_list, σ_list, prob')

# ╔═╡ fb1eeef5-83fc-461d-9881-a7c9584bb1ab
md"### Code 4.19"

# ╔═╡ cb37a0d5-ae6b-4ebc-91a0-144377ef90b7
begin
	indices = collect(Iterators.product(1:length(μ_list), 1:length(σ_list)));
	sample_idx = wsample(vec(indices), vec(prob), 10_000; replace=true)
	sample_μ = μ_list[first.(sample_idx)]
	sample_σ = σ_list[last.(sample_idx)];
end;

# ╔═╡ db6ff059-2a68-4f9e-b592-ce1489569b66
md"### Code 4.20"

# ╔═╡ 708ddb40-ef5d-48c3-952b-1c98431f262f
scatter(sample_μ, sample_σ; alpha=0.1)

# ╔═╡ 6261741f-ff4d-460a-83f2-674553a9c09d
md"### Code 4.21"

# ╔═╡ af99e34c-423b-41da-9fb1-8fda1b1b49c6
let
	p1 = density(sample_μ)
	p2 = density(sample_σ)
	plot(p1, p2, layout=(2,1))
end

# ╔═╡ 1bb443a2-62c9-42de-b508-896e4104794c
md"### Code 4.22"

# ╔═╡ a3fa5e0f-ddd3-4629-b432-174533300e34
hpdi(sample_μ, alpha=0.11)

# ╔═╡ 4779466e-6e1a-41f5-a710-d6de1e6784fd
hpdi(sample_σ, alpha=0.11)

# ╔═╡ 8bacd448-b7d9-4575-9daa-fe4a3ffdd61d
md"### Code 4.23"

# ╔═╡ 8e8284a5-5dfe-4340-9047-c9016d9c8fb8
d3 = sample(howell1_1.height, 20);

# ╔═╡ 4607472d-dd96-473f-ba56-f68a65e25301
md"### Code 4.24"

# ╔═╡ 78f5537e-76e7-44d3-9229-c67850bd615e
begin
	μ_list1 = range(150, 170; length=200)
	σ_list1 = range(4, 20; length=200)

	log_likelihood1 = [
		sum(logpdf.(Normal(μ, σ), d3))
		for μ ∈ μ_list1, σ ∈ σ_list1
	]
	log_prod1 = log_likelihood1 .+ [
		logpdf.(Normal(178, 20), μ) + logpdf.(Uniform(0, 50), σ)
			for μ ∈ μ_list1, σ ∈ σ_list1
	]

	max_prod1 = maximum(log_prod1)
	prob2 = @. exp(log_prod1 - max_prod1)
end

# ╔═╡ b185bfe4-2946-49a2-9d00-27cc5820e31e
begin
	indices1 = collect(Iterators.product(1:length(μ_list1), 1:length(σ_list1)));
	sample2_idx = wsample(vec(indices1), vec(prob2), 10_000; replace=true)
	sample2_μ = μ_list1[first.(sample2_idx)]
	sample2_σ = σ_list1[last.(sample2_idx)]

	scatter(sample2_μ, sample2_σ; alpha=0.1)
end

# ╔═╡ c5758c90-d83b-4860-910f-5c22284fb38a
md"### Code 4.25"

# ╔═╡ 7098d711-dc96-4862-acbb-e0c28ce2287b
let
	density(sample2_σ)
	μ = mean(sample2_σ)
	σ = std(sample2_σ)
	plot!(Normal(μ, σ); label="Normal")
end

# ╔═╡ d950f634-969d-4373-9cf5-313c2ce1e922
md"### Code 4.26"

# ╔═╡ 50040e42-b57c-4587-8c1f-580d83f32ed3
begin
	d = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	d2 = d[d.age .>= 18,:];
end

# ╔═╡ cad3d5d4-075e-4d3b-8e05-e989d07865b7
md"### Code 4.27"

# ╔═╡ 215cf3df-e934-4a45-a1f7-140b12b05c9b
@model function model_height1(height)
    μ ~ Normal(178, 20)
    σ ~ Uniform(0, 50)
    height ~ Normal(μ, σ)
end

# ╔═╡ 6bd4b739-ec92-4585-879a-e086b74ca9bf
md"### Code 4.28"

# ╔═╡ fca9ed55-4bba-41aa-8385-ba11251d8f56
m4_1_2 = sample(model_height1(d2.height), NUTS(), 1000)

# ╔═╡ 1fcd6471-26fd-4ee3-bd60-99a4dc1e6e80
md"### Code 4.29"

# ╔═╡ e93daa5f-621f-4d55-8249-f7f74e0aaad9
describe(m4_1_2; q=[0.055, 0.945])

# ╔═╡ 811aa27c-7ee8-4572-a3af-f94710f42a7f
md"### Code 4.30"

# ╔═╡ 9fa6ae97-566d-4df0-9377-c5316fe43aa5
md"### Code 4.31"

# ╔═╡ 143a02da-ff74-413a-b479-b38223f373b1
@model function m4_2(height)
    μ ~ Normal(178, 0.1)
    σ ~ Uniform(0, 50)
    height ~ Normal(μ, σ)
end

# ╔═╡ 2db2250f-1835-4e5c-a628-441329494648
begin
	init_vals = [mean(d2.height), std(d2.height)]
	chain = sample(m4_2(d2.height), NUTS(), 1000, init_theta=init_vals)
end

# ╔═╡ 7e811d2e-a96b-49f7-a9a6-9ea08d60840e
begin
	m4_2t = m4_2(d2.height)
	post4_2t = sample(m4_2t, NUTS(), 1000)
	describe(post4_2t; q=[0.055, 0.945])
end

# ╔═╡ c3e4e842-a6f4-4ef8-aac6-d7a8308d33c7
md"### Code 4.32"

# ╔═╡ 73de370f-d5e2-405a-adcd-72dfd95e5cab
cov(hcat(m4_1_2[:μ], m4_1_2[:σ]))

# ╔═╡ 546f748d-a82c-4bfe-86d7-fb463de585d7
md"### Code 4.33"

# ╔═╡ 6b065d68-3068-4714-b4f0-46b496541e80
let
	c = cov(hcat(m4_1_2[:μ], m4_1_2[:σ]))
	cov2cor(c, diag(c))
end

# ╔═╡ f57c8744-4fa1-49fd-9f45-857e91ac3762
md"### Code 4.34"

# ╔═╡ 6f2dd2ca-99ee-4464-9db3-5a524679d760
CHNS(m4_1_2)

# ╔═╡ 3945feac-3837-4990-ac87-bf851d8f044b
let
	df = DataFrame(m4_1_2)
	sample(df, 5)
end

# ╔═╡ eca11664-07ac-418c-a0f9-95aa9bc21a5f
md"### Code 4.35"

# ╔═╡ 49103d25-6e41-4236-8ce0-21ac103bb19e
md"### Code 4.36"

# ╔═╡ c94959c6-8058-4014-b630-125405ee09c2
begin
	data = hcat(m4_1_2[:μ], m4_1_2[:σ])
	μ2 = mean(data, dims=1)
	σ2 = cov(data)
	mvn = MvNormal(vec(μ2), σ2)
end

# ╔═╡ 7104075c-e256-499b-8c8b-cc2f50486395
let
	post = rand(mvn, 10_000);
end

# ╔═╡ b48fd6cd-5cf0-4996-a022-f378fae1471b
md" ##### Using Turing's query syntax"

# ╔═╡ 5588d85c-c46a-48ce-9a0a-3d2504b94f90
begin
	@model function m4_3(x, y)
	   s ~ InverseGamma(2, 3)
	   m ~ Normal(0, sqrt(s))
	   x ~ Normal(m, sqrt(s))
	   y ~ Normal(m, sqrt(s))
	end
	
	m4_3t = m4_3(2,4)
	post4_3t = sample(m4_3t, NUTS(), 100)
	
	prob"x=1.0, y=1.0 | model=m4_3t, s = 1.0, m = 1.0"
end

# ╔═╡ 5e652e13-a8ac-4fe5-8024-5ba1a6086219
prob"height=177.0 | model=m4_2t, μ=170.0, σ=20.0"

# ╔═╡ 5c2a7f88-0791-4064-a1ab-ca5effa02c84
prob"x = 1.0, y = 1.0 | chain = post4_3t, model = m4_3t"


# ╔═╡ Cell order:
# ╠═76a7417d-fc80-44c9-b23d-fb75b0e9abba
# ╠═2d3d462e-d6bf-4ebd-ab34-db38f4261a11
# ╠═f09ccff8-4036-4322-90e6-cdceb35eadfe
# ╠═18b5c138-9324-4c5a-86ca-be1e825ed11e
# ╟─a2d950f8-d780-4abe-af37-b189bd4a5e88
# ╠═bbb4a889-77e2-449b-9258-883e7434e641
# ╟─0fb4f0c1-df0a-4734-99c9-9a4671d75bb0
# ╟─5fc53e45-ff63-4542-909a-afd5e2716346
# ╠═e3709463-a936-485b-8173-0f9fda249417
# ╟─58a41999-ae14-423f-b6dd-21fc766374db
# ╠═2f69e9b7-8775-44f9-a493-88d02f8283c7
# ╟─d694159e-7384-465f-b1ab-23ac6b018323
# ╠═08ee5548-ab90-4eba-b151-a37ce6849df6
# ╟─b5f41800-d83a-4f68-80fc-3b52b9b887ea
# ╠═7fa6c286-638d-4a7d-8fef-1dfa9974b67a
# ╟─44b3cfbb-cbb9-4640-84f6-097f52a5b117
# ╠═b70c459c-8bfc-4ece-9e46-e00ca880be89
# ╟─6f8b69ef-acb5-4db4-b2ae-c08f0fe59e6b
# ╟─e9862e8d-047b-4d25-992b-2ed5c9dc26bc
# ╠═4d2231d0-9fa9-4600-b325-ff31da8e327f
# ╟─d34ed551-5cc3-468c-8c17-35f5c92b2280
# ╟─687abeee-b831-4e48-aa8a-9a330bc8d5d1
# ╠═496b257c-ea62-40b2-b1c9-e92fdda1d688
# ╟─4753ae08-1b15-4345-8ae6-a4f876457f93
# ╠═50773db0-910c-42ec-b77a-74ad32f23b55
# ╟─d117f398-6eba-418a-9696-9f34b820f4b5
# ╠═3dd509b8-4e49-42fc-b3bc-5b2aa8c27542
# ╟─48c8de02-bd59-4bca-b103-d4b770b322d7
# ╠═ec60e901-8ea8-4779-bb5d-03d50138e9e6
# ╟─5737278f-8834-4769-9880-1a8ca580f87c
# ╠═90292dea-4a40-47e6-8269-65f76dd1e50a
# ╟─7ce59f8a-3f0d-4562-85dd-911a75d68922
# ╠═fc7fede8-be68-410a-a147-26ce7ab7e008
# ╟─82b5c59d-0277-4212-af32-1f05459e2623
# ╠═ec8fb7f4-83e0-4f3f-8ffa-ec024cfadb88
# ╟─83ae8a58-78af-4044-b476-48bff83f818d
# ╠═a7100c3c-6fab-43b1-af33-a17780f309e6
# ╟─da3e0f8e-9eab-4ebd-bdd7-3738f380fe72
# ╠═dfb5d4fe-d810-4705-b03c-c013b8098315
# ╟─b72c1949-c95f-4a51-91cb-d47c7cedd4cc
# ╠═98ad03fb-0d07-4cf3-9681-2943a8b4e1fd
# ╟─4c904223-defa-4b45-b8b0-e6cf6f9b7cf7
# ╟─8fcc7434-a202-4e42-9473-a66987c3dd18
# ╟─56d515ea-74a6-439a-ab34-be0409a3524e
# ╠═12d0c754-ddc9-4ce5-afba-2b07e57fb826
# ╟─fb1eeef5-83fc-461d-9881-a7c9584bb1ab
# ╠═cb37a0d5-ae6b-4ebc-91a0-144377ef90b7
# ╟─db6ff059-2a68-4f9e-b592-ce1489569b66
# ╠═708ddb40-ef5d-48c3-952b-1c98431f262f
# ╟─6261741f-ff4d-460a-83f2-674553a9c09d
# ╠═af99e34c-423b-41da-9fb1-8fda1b1b49c6
# ╟─1bb443a2-62c9-42de-b508-896e4104794c
# ╠═a3fa5e0f-ddd3-4629-b432-174533300e34
# ╠═4779466e-6e1a-41f5-a710-d6de1e6784fd
# ╟─8bacd448-b7d9-4575-9daa-fe4a3ffdd61d
# ╠═8e8284a5-5dfe-4340-9047-c9016d9c8fb8
# ╟─4607472d-dd96-473f-ba56-f68a65e25301
# ╠═78f5537e-76e7-44d3-9229-c67850bd615e
# ╠═b185bfe4-2946-49a2-9d00-27cc5820e31e
# ╟─c5758c90-d83b-4860-910f-5c22284fb38a
# ╠═7098d711-dc96-4862-acbb-e0c28ce2287b
# ╟─d950f634-969d-4373-9cf5-313c2ce1e922
# ╠═50040e42-b57c-4587-8c1f-580d83f32ed3
# ╟─cad3d5d4-075e-4d3b-8e05-e989d07865b7
# ╠═215cf3df-e934-4a45-a1f7-140b12b05c9b
# ╟─6bd4b739-ec92-4585-879a-e086b74ca9bf
# ╠═fca9ed55-4bba-41aa-8385-ba11251d8f56
# ╟─1fcd6471-26fd-4ee3-bd60-99a4dc1e6e80
# ╠═e93daa5f-621f-4d55-8249-f7f74e0aaad9
# ╟─811aa27c-7ee8-4572-a3af-f94710f42a7f
# ╟─9fa6ae97-566d-4df0-9377-c5316fe43aa5
# ╠═143a02da-ff74-413a-b479-b38223f373b1
# ╠═2db2250f-1835-4e5c-a628-441329494648
# ╠═7e811d2e-a96b-49f7-a9a6-9ea08d60840e
# ╟─c3e4e842-a6f4-4ef8-aac6-d7a8308d33c7
# ╠═73de370f-d5e2-405a-adcd-72dfd95e5cab
# ╟─546f748d-a82c-4bfe-86d7-fb463de585d7
# ╠═6b065d68-3068-4714-b4f0-46b496541e80
# ╟─f57c8744-4fa1-49fd-9f45-857e91ac3762
# ╠═6f2dd2ca-99ee-4464-9db3-5a524679d760
# ╠═3945feac-3837-4990-ac87-bf851d8f044b
# ╟─eca11664-07ac-418c-a0f9-95aa9bc21a5f
# ╟─49103d25-6e41-4236-8ce0-21ac103bb19e
# ╠═c94959c6-8058-4014-b630-125405ee09c2
# ╠═7104075c-e256-499b-8c8b-cc2f50486395
# ╟─b48fd6cd-5cf0-4996-a022-f378fae1471b
# ╠═5588d85c-c46a-48ce-9a0a-3d2504b94f90
# ╠═5e652e13-a8ac-4fe5-8024-5ba1a6086219
# ╠═5c2a7f88-0791-4064-a1ab-ca5effa02c84
