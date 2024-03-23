### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ bddf36ac-5547-491c-94bb-763ce7df3cd9
using Pkg, DrWatson

# ╔═╡ 7bf4154b-0009-46e4-9d03-8e3b81a202f2
begin
	using Optim
	using CSV
	using Random
	using StatsBase
	using DataFrames
	using Turing
	using StatsPlots
	using StatsFuns
	using LaTeXStrings
	using StatisticalRethinking
	using StatisticalRethinking: link
	using StatisticalRethinkingPlots
	using ParetoSmooth
	using ParetoSmoothedImportanceSampling
	using Logging
end

# ╔═╡ 8720cc1f-58a0-4a70-9497-480456b26f85
begin
	default(labels=false)
	Logging.disable_logging(Logging.Warn);
end

# ╔═╡ fdbf8fe2-7e99-4da6-bd8d-6a21bb68bc71
md"## 8.1 Building an interaction."

# ╔═╡ 04affaff-2593-4b6c-8638-ab54fde6d2cc
md"### Code 8.1"

# ╔═╡ 73ee2643-330f-485f-94f5-be46fbd33d4f
begin
	rugged = CSV.read(sr_datadir("rugged.csv"), DataFrame)
	dd = rugged[completecases(rugged, :rgdppc_2000),:]
	dd[:,:log_gdp] = log.(dd.rgdppc_2000);
	dd[:,:log_gdp_std] = dd.log_gdp / mean(dd.log_gdp)
	dd[:,:rugged_std] = dd.rugged / maximum(dd.rugged)
end;

# ╔═╡ eb3fef03-acf9-4c41-9760-e0f492df55e2
md"### Code 8.2"

# ╔═╡ 3eeab20c-6d42-4b31-bd8f-bfa691cb4587
begin
	r̄ = mean(dd.rugged_std)
	
	@model function model_m8_1(rugged_std, log_gdp_std)
	    σ ~ Exponential()
	    a ~ Normal(1, 1)
	    b ~ Normal(0, 1)
	    μ = @. a + b * (rugged_std - r̄)
	    log_gdp_std ~ MvNormal(μ, σ)
	end
end

# ╔═╡ bc0a3a04-1b8c-4332-a6ba-3da2b8eaea05
md"### Code 8.3"

# ╔═╡ f74d5bf7-ae57-4678-be5f-25af1959173b
begin
	m8_1_p = sample(model_m8_1(dd.rugged_std, dd.log_gdp_std), Prior(), 1000)
	m8_1_p_df = DataFrame(m8_1_p);
end

# ╔═╡ 1bc733fa-0062-4521-8617-4e757a9fff24
begin
	rugged_seq = range(-0.1, 1.1; length=30)
	μ = link(m8_1_p_df, (r, x) -> r.a + r.b*(x - r̄), rugged_seq)
	μ = hcat(μ...)
	
	p = plot(
	    xlim=(0, 1),
	    ylim=(0.5, 1.5),
	    #title=L"a \sim \mathcal{N}(1,1), b \sim \mathcal{N}(0, 1)",
	    xlab="ruggedness", ylab="log GDP",
	)
	hline!(collect(extrema(dd.log_gdp_std)); c=:black, s=:dash)
	for μ₀ ∈ first(eachrow(μ), 50)
	    plot!(rugged_seq, μ₀; c=:black, alpha=0.3)
	end
	p
end

# ╔═╡ 4f444d5f-6058-4cb8-975f-17782d8b86c9
md"### Code 8.4"

# ╔═╡ cee469b4-83b1-4e29-aa75-8bb6dfc4d207
mean(abs.(m8_1_p_df.b) .> 0.6)

# ╔═╡ 3d677fac-ff09-45b1-9377-4d29cfa86a7f
md"### Code 8.5"

# ╔═╡ 8f904830-dc3e-4875-9af7-1e7285e7fc43
r̄_std = mean(dd.rugged_std)

# ╔═╡ ad981e3c-df35-43f0-9b19-e128a9f2a12d
@model function model_m8_1a(rugged_std, log_gdp_std)
    σ ~ Exponential()
    a ~ Normal(1, 0.1)
    b ~ Normal(0, 0.3)
    μ = @. a + b * (rugged_std - r̄_std)
    log_gdp_std ~ MvNormal(μ, σ)
end

# ╔═╡ 527fd7d7-7375-4dee-9405-cee0363ca92e
begin
	m8_1 = sample(model_m8_1a(dd.rugged_std, dd.log_gdp_std), NUTS(), 1000)
	m8_1_df = DataFrame(m8_1)
end

# ╔═╡ d11c957d-bc78-4011-9dc8-5335a9d3f4f7
md"### Code 8.6"

# ╔═╡ 16574ab6-a53e-4f9d-bbbe-23714ee39abd
describe(m8_1_df)

# ╔═╡ 40799cc2-cdaa-4a86-a678-1b50b3dff6cc
md"### Code 8.7"

# ╔═╡ ca390706-13d7-4aa0-952e-8f06b4bc370e
dd[:,:cid] = @. ifelse(dd.cont_africa == 1, 1, 2);

# ╔═╡ 1d902b61-2387-4664-a39e-7ff73440dd43
md"### Code 8.8"

# ╔═╡ 62e35bf2-cefd-413d-a2ea-60e7193c94b6
@model function model_m8_2(rugged_std, cid,  log_gdp_std)
    σ ~ Exponential()
    a ~ MvNormal([1, 1], 0.1)
    b ~ Normal(0, 0.3)
    μ = @. a[cid] + b * (rugged_std - r̄)
    log_gdp_std ~ MvNormal(μ, σ)
end

# ╔═╡ 43e22387-2afd-4ed1-80c2-608f82c224f6
begin
	m8_2 = sample(model_m8_2(dd.rugged_std, dd.cid,  dd.log_gdp_std), NUTS(), 1000)
	m8_2_df = DataFrame(m8_2);
end# -

# ╔═╡ f5cb4db3-99d4-4498-ba6d-5b36064b5328
md"### Code 8.9"

# ╔═╡ 5f24eb52-4bf3-47b9-ab55-f92aaeca0d45
let
	# Compute log likelihoods for both models
	fun = (r, (x,y)) -> normlogpdf(r.a + r.b * (x - r̄), r.σ, y)
	global m8_1_ll = link(m8_1_df, fun, zip(dd.rugged_std, dd.log_gdp_std))
	m8_1_ll = hcat(m8_1_ll...)
end

# ╔═╡ 612456bd-2bee-4f1c-beed-d3f697c80509
let
	# need DF with a as a vector of both a[1] and a[2]
	global df = DataFrame(m8_2_df)
	df[!,:a] = collect.(zip(m8_2_df.:"a[1]", m8_2_df.:"a[2]"))
	
	fun = (r, (x,c,y)) -> normlogpdf(r.a[c] + r.b * (x - r̄), r.σ, y)
	global m8_2_ll = link(df, fun, zip(dd.rugged_std, dd.cid, dd.log_gdp_std))
	m8_2_ll = hcat(m8_2_ll...);
	
	compare([m8_1_ll, m8_2_ll], :waic, mnames=["m8.1", "m8.2"])
end

# ╔═╡ 8ca7b83a-ffab-42cc-bca1-13cb60a45b97
md"### Code 8.10"

# ╔═╡ fd2ed8f7-c276-4df6-8d08-9d4ea02effd7
describe(m8_2_df)

# ╔═╡ e396ec23-f67f-4f1b-b6e1-3968120463a9
md"### Code 8.11"

# ╔═╡ 38d6aeb9-5b1f-456b-9de6-7484ba94f9c8
PI(map(r -> r[1] - r[2], df.a))

# ╔═╡ 4850e051-166a-4b4a-b7bd-4181578338dc
md"### Code 8.12"

# ╔═╡ 356ef2ab-ee33-42ba-802e-84e1d7425a61
begin
	rugged_seq_2 = range(-0.1, 1.1, length=30)
	africa     = link(df, (r, x) -> r.a[1] + r.b*(x-r̄), rugged_seq_2)
	africa     = hcat(africa...)'
	not_africa = link(df, (r, x) -> r.a[2] + r.b*(x-r̄), rugged_seq_2)
	not_africa = hcat(not_africa...)'
	
	μₐ = mean.(eachrow(africa))
	μₙ = mean.(eachrow(not_africa))
	PIₐ = PI.(eachrow(africa))
	PIₐ = vcat(PIₐ'...)
	PIₙ = PI.(eachrow(not_africa))
	PIₙ = vcat(PIₙ'...);
end

# ╔═╡ 02a727d7-5095-4728-abd4-413100d56acd
let
	p = plot(xlab="ruggedness (std)", ylab="log GDP")
	scatter!(dd.rugged_std[dd.cid.==1], dd.log_gdp_std[dd.cid.==1], c=:blue)
	scatter!(dd.rugged_std[dd.cid.==2], dd.log_gdp_std[dd.cid.==2], c=:white)
	
	plot!(rugged_seq, [μₐ, μₐ], c=:black, fillrange=PIₐ, fillalpha=0.4)
	plot!(rugged_seq, [μₙ, μₙ], c=:black, fillrange=PIₙ, fillalpha=0.2)
	annotate!([
	    (1, 0.87, ("Africa", 9)),
	    (1, 1.05, ("Not Africa", 9))
	])
end

# ╔═╡ 083b4a77-da6e-43e7-a7f0-f45cb2b56ff8
md"### Code 8.13"

# ╔═╡ 3571a823-63f5-4124-ab8f-d7b060b82d6a
@model function model_m8_3(rugged_std, cid,  log_gdp_std)
    σ ~ Exponential()
    a ~ MvNormal([1, 1], 0.1)
    b ~ MvNormal([0, 0], 0.3)
    μ = @. a[cid] + b[cid] * (rugged_std - r̄)
    log_gdp_std ~ MvNormal(μ, σ)
end

# ╔═╡ f2a34af7-2706-4d58-a3c9-88ca3539538e
begin
	m8_3 = sample(model_m8_3(dd.rugged_std, dd.cid,  dd.log_gdp_std), NUTS(), 1000)
	m8_3_df = DataFrame(m8_3);
end

# ╔═╡ e77013bc-1386-4f07-bf64-2571fdbc5212
md"### Code 8.14"

# ╔═╡ 87311d6b-db5f-41d3-a2ec-1e76aa13dd8d
describe(m8_3_df)

# ╔═╡ 1462e546-c794-45aa-8ac2-097b5a0fc168
md"### Code 8.15"

# ╔═╡ 3a55d73d-2cf8-4411-a075-81157a544973
let
	global df3 = DataFrame(m8_3_df)
	df3[!,:a] = collect.(zip(m8_3_df.:"a[1]", m8_3_df.:"a[2]"))
	df3[!,:b] = collect.(zip(m8_3_df.:"b[1]", m8_3_df.:"b[2]"))
	
	fun = (r, (x,c,y)) -> normlogpdf(r.a[c] + r.b[c] * (x - r̄), r.σ, y)
	global m8_3_ll = link(df3, fun, zip(dd.rugged_std, dd.cid, dd.log_gdp_std))
	m8_3_ll = hcat(m8_3_ll...);
	
	compare([m8_1_ll, m8_2_ll, m8_3_ll], :psis, mnames=["m8.1", "m8.2", "m8.3"])
end

# ╔═╡ 80eca5f8-762b-4b73-a246-5883a04d72e5
md"### Code 8.16"

# ╔═╡ 9b1a7350-4311-48b2-91ce-8f4f87e34aef
let
	t = m8_3_ll'
	m8_3_t = collect(reshape(t, size(t)..., 1))
	PSIS_m8_3 = psis_loo(m8_3_t)
	scatter(PSIS_m8_3.pointwise(:pareto_k))
end

# ╔═╡ b384f828-985d-4f6e-9040-2af4bc888ec8
md"### Code 8.17"

# ╔═╡ 8d2cc03b-4e3f-409d-8bdb-e99ed0365046
let
	# build data
	africa     = link(df3, (r, x) -> r.a[1] + r.b[1]*(x-r̄), rugged_seq)
	africa     = hcat(africa...)'
	not_africa = link(df3, (r, x) -> r.a[2] + r.b[2]*(x-r̄), rugged_seq)
	not_africa = hcat(not_africa...)'
	
	μₐ = mean.(eachrow(africa))
	μₙ = mean.(eachrow(not_africa))
	PIₐ = PI.(eachrow(africa))
	PIₐ = vcat(PIₐ'...)
	PIₙ = PI.(eachrow(not_africa))
	PIₙ = vcat(PIₙ'...);
	
	# plot Africa, cid=1
	p1 = plot(xlab="ruggedness (std)", ylab="log GDP", title="African nations")
	scatter!(dd.rugged_std[dd.cid.==1], dd.log_gdp_std[dd.cid.==1], c=:blue)
	plot!(rugged_seq, [μₐ, μₐ], c=:blue, fillrange=PIₐ, fillalpha=0.2)
	
	# plot non Africa, cid=2
	p2 = plot(xlab="ruggedness (std)", ylab="log GDP", title="Non-African nations")
	scatter!(dd.rugged_std[dd.cid.==2], dd.log_gdp_std[dd.cid.==2], c=:white)
	plot!(rugged_seq, [μₙ, μₙ], c=:black, fillrange=PIₙ, fillalpha=0.2)
	
	plot(p1, p2, size=(800, 400))
end

# ╔═╡ e8cbf4c4-885e-4e2f-b058-4b41a3ce6857
md"## 8.2 Symmetry of interations"

# ╔═╡ b0a849d8-315a-44ea-ac6a-769bf47de6ff
md"### Code 8.18"

# ╔═╡ 3383537c-7af4-4de9-b9eb-9389a02d4aaf
let
	rugged_seq = range(-0.2, 1.2, length=30)
	μA = link(df3, (r, x) -> r.a[1] + r.b[1]*(x-r̄), rugged_seq)
	μA = vcat(μA'...)
	μN = link(df3, (r, x) -> r.a[2] + r.b[2]*(x-r̄), rugged_seq)
	μN = vcat(μN'...)
	delta = μA .- μN;
	
	# +
	μ = mean.(eachrow(delta))
	PI_v = PI.(eachrow(delta))
	PI_v = vcat(PI_v'...)
	
	plot(xlab="ruggedness", ylab="expected difference log GDP",)
	plot!(rugged_seq, [μ, μ], c=:blue, fillrange=PI_v, fillalpha=0.2)
	hline!([0.0], s=:dash, c=:black)
	annotate!([
	    (0.0, 0.03, ("Africa higher GPD", 10)),
	    (0.0, -0.03, ("Africa lower GPD", 10)),
	])
end

# ╔═╡ 7d5e4157-7abe-4fb1-893a-51a520a934d0
md"## 8.3 Continuous interaction"

# ╔═╡ 48573194-9cf9-4c9d-b883-62547e4debbc
md"### Code 8.19"

# ╔═╡ 9025cb4c-80cc-463f-b588-913e59c35560
begin
	tulips = CSV.read(sr_datadir("tulips.csv"), DataFrame)
	describe(tulips)
end

# ╔═╡ b2dbe2bd-19ae-4a8f-a408-0353cccd3b84
md"## Code 8.20"

# ╔═╡ 9c905e24-fa5d-4169-bf73-684f5e0a5b8c
begin
	tulips.blooms_std = tulips.blooms / maximum(tulips.blooms)
	tulips.water_cent = tulips.water .- mean(tulips.water)
	tulips.shade_cent = tulips.shade .- mean(tulips.shade);
end;

# ╔═╡ 2f086d7b-658d-42e8-a112-d22807a66586
md"### Code 8.21"

# ╔═╡ 69e42d8d-062a-4061-896b-2781a4aab5b0
let
	Random.seed!(1)
	a = rand(Normal(0.5, 1), 10^4)
	sum(@. (a < 0) | (a > 1))/length(a)
end

# ╔═╡ 219bc61b-a458-488c-91cb-0b3708c74a48
md"### Code 8.22"

# ╔═╡ 4443c25d-7f6a-4d01-b672-1c9606ab28f8
let
	Random.seed!(1)
	a = rand(Normal(0.5, 0.25), 10^4)
	sum(@. (a < 0) | (a > 1))/length(a)
end

# ╔═╡ 738c7fb8-8f9c-4fce-ad6a-b58fd7085be6
md"### Code 8.23"

# ╔═╡ 01b45f6a-adbe-4791-88d5-fa0623bb107c
# +
@model function m8_4(water_cent, shade_cent, blooms_std)
    a ~ Normal(0.5, 0.25)
    bw ~ Normal(0, 0.25)
    bs ~ Normal(0, 0.25)
    μ = @. a + bw*water_cent + bs*shade_cent
    σ ~ Exponential(1)
    blooms_std ~ MvNormal(μ, σ)
end

# ╔═╡ 715b797a-88c6-46b9-87b9-a59d5765cd99
begin
	m8_4_c = sample(m8_4(tulips.water_cent, tulips.shade_cent, tulips.blooms_std), 
		NUTS(), 1000)
	m8_4_df = DataFrame(m8_4_c)
	describe(m8_4_df)
end

# ╔═╡ 1f7c87b9-a93c-442e-b5b2-8501a289200c
md"### Code 8.24"

# ╔═╡ 6a5a4d0f-5932-4862-87b0-1fafed29318e
@model function m8_5(water_cent, shade_cent, blooms_std)
    a ~ Normal(0.5, 0.25)
    bw ~ Normal(0, 0.25)
    bs ~ Normal(0, 0.25)
    bws ~ Normal(0, 0.25)
    μ = @. a + bw*water_cent + bs*shade_cent + bws*water_cent*shade_cent
    σ ~ Exponential(1)
    blooms_std ~ MvNormal(μ, σ)
end

# ╔═╡ a369c511-017d-406f-9178-b4cced1b97e9
begin
	m8_5_c = sample(m8_5(tulips.water_cent, tulips.shade_cent, tulips.blooms_std),
		NUTS(), 1000)
	m8_5_df = DataFrame(m8_5_c)
	describe(m8_5_df)
end

# ╔═╡ 04be14ef-fe26-447b-823e-243d9cba8466
md"### Code 8.25"

# ╔═╡ f822f689-a92c-4535-ae88-ef604e31753f
let
	plts = []
	
	for shade ∈ -1:1
	    idx = findall(==(shade), tulips.shade_cent)
	    p = plot(xlims=(-1.2,1.2), ylims=(-.2,1.2), xlab="water", ylab="blooms", 
	             title="shade=$shade", titlefontsize=12)
	    scatter!(tulips.water_cent[idx], tulips.blooms_std[idx])
	    water_seq = -1:1
	    mu = link(m8_4_df, (r, water) -> r.a + r.bw * water + r.bs * shade, 
			water_seq)
	    mu = hcat(mu...);
	    for μ ∈ first(eachrow(mu), 20)
	        plot!(water_seq, μ, c=:black, alpha=0.2)
	    end
	    push!(plts, p)
	end
	plot(plts..., layout=(1, 3), size=(800, 400), plot_title="m8.4 post", 
		plot_titlefontsize=14)
end

# ╔═╡ 745a0c84-4243-4efa-9871-7d2dcd096553
let
	plts = []
	
	for shade ∈ -1:1
	    idx = findall(==(shade), tulips.shade_cent)
	    p = plot(xlims=(-1.2,1.2), ylims=(-.2,1.2), xlab="water", ylab="blooms", 
	             title="shade=$shade", titlefontsize=12)
	    scatter!(tulips.water_cent[idx], tulips.blooms_std[idx])
	    water_seq = -1:1
	    mu = link(m8_5_df, (r, water) -> r.a + r.bw*water + r.bs*shade + 
			r.bws*water*shade, water_seq)
	    mu = hcat(mu...);
	    for μ ∈ first(eachrow(mu), 20)
	        plot!(water_seq, μ, c=:black, alpha=0.2)
	    end
	    push!(plts, p)
	end
	plot(plts..., layout=(1, 3), size=(800, 400), plot_title="m8.5 post", 
		plot_titlefontsize=14)
end

# ╔═╡ e0bc4345-c125-42ca-a088-b70a7c2f762c
md"### Code 8.26"

# ╔═╡ 0c58026c-ec86-4b16-9f00-398bed8d627a
begin
	Random.seed!(7)
	m8_5p_c = sample(m8_5(tulips.water_cent, tulips.shade_cent, tulips.blooms_std), 
		Prior(), 1000)
	m8_5p_df = DataFrame(m8_5p_c);
end

# ╔═╡ 20a3b233-a892-4518-8912-a49dc8ee2126
let
	plts = []
	
	for shade ∈ -1:1
	    p = plot(xlims=(-1, 1), ylims=(-0.5, 1.5), xlab="water", ylab="blooms", 
	             title="shade=$shade", titlefontsize=12)
	    water_seq = -1:1
	    mu = link(m8_5p_df, (r, water) -> r.a + r.bw*water + r.bs*shade + 
			r.bws*water*shade, water_seq)
	    mu = hcat(mu...);
	    for μ ∈ first(eachrow(mu), 20)
	        plot!(water_seq, μ, c=:black, alpha=0.2)
	    end
	    hline!([0.0, 1.0], s=:dash, c=:black)
	    push!(plts, p)
	end
	plot(plts..., layout=(1, 3), size=(800, 400), plot_title="m8.5 prior", 
		plot_titlefontsize=14)
end

# ╔═╡ Cell order:
# ╠═bddf36ac-5547-491c-94bb-763ce7df3cd9
# ╠═7bf4154b-0009-46e4-9d03-8e3b81a202f2
# ╠═8720cc1f-58a0-4a70-9497-480456b26f85
# ╟─fdbf8fe2-7e99-4da6-bd8d-6a21bb68bc71
# ╟─04affaff-2593-4b6c-8638-ab54fde6d2cc
# ╠═73ee2643-330f-485f-94f5-be46fbd33d4f
# ╟─eb3fef03-acf9-4c41-9760-e0f492df55e2
# ╠═3eeab20c-6d42-4b31-bd8f-bfa691cb4587
# ╟─bc0a3a04-1b8c-4332-a6ba-3da2b8eaea05
# ╠═f74d5bf7-ae57-4678-be5f-25af1959173b
# ╠═1bc733fa-0062-4521-8617-4e757a9fff24
# ╟─4f444d5f-6058-4cb8-975f-17782d8b86c9
# ╠═cee469b4-83b1-4e29-aa75-8bb6dfc4d207
# ╟─3d677fac-ff09-45b1-9377-4d29cfa86a7f
# ╠═8f904830-dc3e-4875-9af7-1e7285e7fc43
# ╠═ad981e3c-df35-43f0-9b19-e128a9f2a12d
# ╠═527fd7d7-7375-4dee-9405-cee0363ca92e
# ╠═d11c957d-bc78-4011-9dc8-5335a9d3f4f7
# ╠═16574ab6-a53e-4f9d-bbbe-23714ee39abd
# ╟─40799cc2-cdaa-4a86-a678-1b50b3dff6cc
# ╠═ca390706-13d7-4aa0-952e-8f06b4bc370e
# ╠═1d902b61-2387-4664-a39e-7ff73440dd43
# ╠═62e35bf2-cefd-413d-a2ea-60e7193c94b6
# ╠═43e22387-2afd-4ed1-80c2-608f82c224f6
# ╟─f5cb4db3-99d4-4498-ba6d-5b36064b5328
# ╠═5f24eb52-4bf3-47b9-ab55-f92aaeca0d45
# ╠═612456bd-2bee-4f1c-beed-d3f697c80509
# ╟─8ca7b83a-ffab-42cc-bca1-13cb60a45b97
# ╠═fd2ed8f7-c276-4df6-8d08-9d4ea02effd7
# ╟─e396ec23-f67f-4f1b-b6e1-3968120463a9
# ╠═38d6aeb9-5b1f-456b-9de6-7484ba94f9c8
# ╟─4850e051-166a-4b4a-b7bd-4181578338dc
# ╠═356ef2ab-ee33-42ba-802e-84e1d7425a61
# ╠═02a727d7-5095-4728-abd4-413100d56acd
# ╠═083b4a77-da6e-43e7-a7f0-f45cb2b56ff8
# ╠═3571a823-63f5-4124-ab8f-d7b060b82d6a
# ╠═f2a34af7-2706-4d58-a3c9-88ca3539538e
# ╠═e77013bc-1386-4f07-bf64-2571fdbc5212
# ╠═87311d6b-db5f-41d3-a2ec-1e76aa13dd8d
# ╟─1462e546-c794-45aa-8ac2-097b5a0fc168
# ╠═3a55d73d-2cf8-4411-a075-81157a544973
# ╟─80eca5f8-762b-4b73-a246-5883a04d72e5
# ╠═9b1a7350-4311-48b2-91ce-8f4f87e34aef
# ╠═b384f828-985d-4f6e-9040-2af4bc888ec8
# ╠═8d2cc03b-4e3f-409d-8bdb-e99ed0365046
# ╠═e8cbf4c4-885e-4e2f-b058-4b41a3ce6857
# ╠═b0a849d8-315a-44ea-ac6a-769bf47de6ff
# ╠═3383537c-7af4-4de9-b9eb-9389a02d4aaf
# ╟─7d5e4157-7abe-4fb1-893a-51a520a934d0
# ╟─48573194-9cf9-4c9d-b883-62547e4debbc
# ╠═9025cb4c-80cc-463f-b588-913e59c35560
# ╟─b2dbe2bd-19ae-4a8f-a408-0353cccd3b84
# ╠═9c905e24-fa5d-4169-bf73-684f5e0a5b8c
# ╟─2f086d7b-658d-42e8-a112-d22807a66586
# ╠═69e42d8d-062a-4061-896b-2781a4aab5b0
# ╟─219bc61b-a458-488c-91cb-0b3708c74a48
# ╠═4443c25d-7f6a-4d01-b672-1c9606ab28f8
# ╟─738c7fb8-8f9c-4fce-ad6a-b58fd7085be6
# ╠═01b45f6a-adbe-4791-88d5-fa0623bb107c
# ╠═715b797a-88c6-46b9-87b9-a59d5765cd99
# ╟─1f7c87b9-a93c-442e-b5b2-8501a289200c
# ╠═6a5a4d0f-5932-4862-87b0-1fafed29318e
# ╠═a369c511-017d-406f-9178-b4cced1b97e9
# ╟─04be14ef-fe26-447b-823e-243d9cba8466
# ╠═f822f689-a92c-4535-ae88-ef604e31753f
# ╠═745a0c84-4243-4efa-9871-7d2dcd096553
# ╟─e0bc4345-c125-42ca-a088-b70a7c2f762c
# ╠═0c58026c-ec86-4b16-9f00-398bed8d627a
# ╠═20a3b233-a892-4518-8912-a49dc8ee2126
