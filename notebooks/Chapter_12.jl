### A Pluto.jl notebook ###
# v0.19.3

using Markdown
using InteractiveUtils

# ╔═╡ 2b447c25-cc6e-4b40-a18a-fcab8eb64d8e
using Pkg, DrWatson

# ╔═╡ c9ea9876-b966-41b4-ac02-a6b1040e7119
begin
	using Turing
	using DataFrames
	using CSV
	using Random
	using Distributions
	using StatisticalRethinking
	using StatisticalRethinking: link
	using StatisticalRethinkingPlots
	using ParetoSmooth
	using StatsPlots
	using StatsBase
	using FreqTables
	using Logging
end

# ╔═╡ 147708bd-6ded-4316-aca0-1f6524ba3c53
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

# ╔═╡ 2c26b000-36ee-473a-a7c5-d6a746096ec9
begin
	default(label=false);
	Logging.disable_logging(Logging.Warn);
end

# ╔═╡ 7c87ccc8-8b76-45f1-b10f-b04bba307ce9
md" ## 12.1 Over-dispersed counts"

# ╔═╡ 5be273e4-faae-4f06-934c-a1da019e8b89
md" #### Code 12.1"

# ╔═╡ 1651ab99-e88a-40b7-9092-daf7f31e5604
begin
	# define alias for Beta(α, β), see: https://en.wikipedia.org/wiki/Beta_distribution#Mean_and_sample_size
	Beta2(μ, ν) = Beta(μ*ν, (1-μ)*ν)
	BetaBinomial2(n, μ, ν) = BetaBinomial(n, μ*ν, (1-μ)*ν)
end

# ╔═╡ fdd7dad1-facd-4ce3-9a72-20ce2c4cbe6a
let
	p̄ = 0.5
	θ = 5
	plot(Beta2(p̄, θ), xlab="probability", ylab="Density")
end

# ╔═╡ 86d02852-70e2-4fcc-bc8c-8d2b881fa248
md" #### Code 12.2"

# ╔═╡ f2af75cc-361e-4a45-a1b4-14a4ce404fa3
begin
	ucbadmit = CSV.read(sr_datadir("UCBadmit.csv"), DataFrame)
	ucbadmit.gid = @. ifelse(ucbadmit.gender == "male", 1, 2)
end;

# ╔═╡ 8d1b576f-e2a4-423d-80cc-0e69f6a376d4
@model function m12_1(A, N, gid)
    a ~ MvNormal([0, 0], 1.5)
    p̄ = @. logistic(a[gid])
    ϕ ~ Exponential(1)
    θ = ϕ + 2
    @. A ~ BetaBinomial2(N, p̄, θ)
end

# ╔═╡ 9ce82e0e-bdcf-426c-b87e-32346c658b87
begin
	Random.seed!(1)
	m12_1_ch = sample(m12_1(ucbadmit.admit, ucbadmit.applications, ucbadmit.gid), NUTS(), 1000)
	m12_1_df = DataFrame(m12_1_ch);
end

# ╔═╡ b339874c-e364-4f8c-bd62-5d66b175f0be
md" #### Code 12.3"

# ╔═╡ 50661b7a-3484-4c52-b87c-4f44f119e164
let
	m12_1_df.θ = m12_1_df.ϕ .+ 2
	m12_1_df.da = m12_1_df."a[1]" .- m12_1_df."a[2]"
	PRECIS(m12_1_df)
end

# ╔═╡ 4335575b-2ce8-41f0-a7ed-295789a5b723
md" #### Code 12.4"

# ╔═╡ c8959699-8886-423a-b994-3180c639de1b
let
	gid = 2
	p̄ = m12_1_df[:, "a[$gid]"] .|> logistic |> mean
	θ = mean(m12_1_df.θ)
	plot(Beta2(p̄, θ), lw=2, c=:black, ylab="Density", 
	    xlab="probability admit", ylims=(0, 3))
	
	for (a, θ) ∈ first(zip(m12_1_df[:, "a[$gid]"], m12_1_df.θ), 50)
	    plot!(Beta2(logistic(a), θ), c=:black, alpha=0.2)
	end
	gender = gid == 2 ? "female" : "male"
	title!("distribution of $gender admission rates")
end

# ╔═╡ 7959fa22-560e-44d0-96ee-d2e7bf51cbf9
md" #### Code 12.5"

# ╔═╡ 7a9ad8dc-019d-4a32-b8b3-ee5a2703e609
md" #### Code 12.6"

# ╔═╡ 5e3a0ec2-5066-4ab6-881d-af6113d23974
begin
	kline = CSV.read(sr_datadir("Kline.csv"), DataFrame)
	kline.P = standardize(ZScoreTransform, log.(kline.population))
	kline.contact_id = ifelse.(kline.contact .== "high", 2, 1)
end;

# ╔═╡ 29167add-b592-4a1e-8742-88fa3264286d
@model function m12_2(T, P, cid)
    g ~ Exponential()
    ϕ ~ Exponential()
    a ~ MvNormal([1,1], 1)
    b₁ ~ Exponential()
    b₂ ~ Exponential()
    b = [b₁, b₂]
    λ = @. exp(a[cid])*(P^b[cid]) / g
    p = 1/(ϕ+1)
    r = λ/ϕ
    clamp!(r, 0.01, Inf)
    p = clamp(p, 0.01, 1)
    @. T ~ NegativeBinomial(r, p)
end

# ╔═╡ cc4b5a4e-c2c6-4295-b4ff-8f73e84fcd88
begin
	Random.seed!(1)
	m12_2_ch = sample(m12_2(kline.total_tools, kline.population, kline.contact_id), NUTS(), 1000)
	m12_2_df = DataFrame(m12_2_ch)
	PRECIS(m12_2_df)
end

# ╔═╡ 512a3b1f-2e79-4777-b976-6d310725e9f8
md" ## 12.2 Zero-inflated outcomes"

# ╔═╡ 91de0472-61b8-47f4-babb-e07539fa3a55
md" #### Codes 12.7 & 12.8"

# ╔═╡ dc4292b2-5aa8-4413-b4a3-91456f1d1d2a
md" #### Code 12.9"

# ╔═╡ 3288b258-640b-414f-90d4-0bc1517aa506
# Based on this discussion
# https://github.com/StatisticalRethinkingJulia/SR2TuringPluto.jl/issues/1

# ╔═╡ 45e9dc3c-bd28-4eaf-a56a-d12f4804f341
struct ZIPoisson{T1,T2} <: DiscreteUnivariateDistribution
    λ::T1
    w::T2
end

# ╔═╡ ab414dc7-7f31-4b0a-9d65-1268940fbc50
function logpdf(d::ZIPoisson, y::Int)
    if y == 0
        logsumexp([log(d.w), log(1 - d.w) - d.λ])
    else
        log(1 - d.w) + logpdf(Poisson(d.λ), y)
    end
end

# ╔═╡ c27a8622-edca-4f19-91d0-1682ebc9cc7d
fun = (r, (N, gid)) -> begin
    p̄ = logistic(get(r, "a[$gid]", 0))
    rand(BetaBinomial2(N, p̄, r.θ))
end

# ╔═╡ c799de28-efa8-4466-812f-ac3c79b5d835
let
	#import Distributions: logpdf, rand

	Random.seed!(1)
	adm_rate = ucbadmit.admit ./ ucbadmit.applications
	pred_adm = link(m12_1_df, fun, zip(ucbadmit.applications, ucbadmit.gid))
	pred_rates = pred_adm ./ ucbadmit.applications

	μ_adm = mean.(pred_rates)
	σ = std.(pred_rates) ./ 2
	ci_adm = PI.(pred_rates)

	scatter(adm_rate, xlab="case", ylab="A", title="Posterior validation check")
	scatter!(μ_adm, mc=:white, yerr=σ)
	scatter!(first.(ci_adm), shape=:cross, c=:black)
	scatter!(last.(ci_adm), shape=:cross, c=:black)
end

# ╔═╡ 15d28187-c155-4332-bde4-90c2d40b7a36
let
	# define parameters
	prob_drink = 0.2 # 20% of days
	rate_work = 1    # average 1 manuscript per day

	# sample one year of production
	N = 365

	# simulate days monks drink
	Random.seed!(365)
	drink = rand(Binomial(1, prob_drink), N)

	# simulate manuscripts completed
	y = (1 .- drink).*rand(Poisson(rate_work), N);
	p = histogram(y, xlab="manuscripts completed", ylab="Frequency")
	zeros_drink = sum(drink)
	zeros_work = sum(@. (y == 0) & (drink == 0))
	bar!([0], [zeros_work], bar_width=0.3)
	p
end

# ╔═╡ 73d535fd-227a-4c4b-a0ca-ed926ecf86c6
@model function m12_3(y)
    ap ~ Normal(-1.5, 1)
    al ~ Normal(1, 0.5)
    λ = exp(al)
    p = logistic(ap)
    y .~ ZIPoisson(λ, p)
end

# ╔═╡ cb63530d-4041-4a53-92c7-73a5c5705564
begin
	m12_3_ch = sample(m12_3(y), NUTS(), 1000)
	m12_3_df = DataFrame(m12_3_ch)
	precis(m12_3_df)
end

# ╔═╡ 3b4f512c-7fc3-4eda-90a3-57d7ed313f87
md" #### Code 12.10"

# ╔═╡ 719b6fe8-03e0-4e1d-8be8-6d642f78e86d
let
	m12_3_df.ap .|> logistic |> mean,
	exp.(m12_3_df.al) |> mean
end

# ╔═╡ e3ea5e44-b189-4c2e-b7d0-5bdc26477ef8
md" ## 12.3 Ordered categorical outcomes"

# ╔═╡ 8c89b94a-3250-4abf-8eda-35fd2b42f708
md" #### Code 12.12"

# ╔═╡ a6bcc9c6-4be7-4390-9fa3-6a124ed8afc1
begin
	trolley = CSV.read(sr_datadir("Trolley.csv"), DataFrame)
	describe(trolley)
end

# ╔═╡ cbf2db2b-5e72-4ad1-9fc3-0944c17db2df
md" #### Code 12.13"

# ╔═╡ c1943565-d495-480d-86ce-ea6f4e760b53
histogram(trolley.response, xlab="response")

# ╔═╡ ee421d11-e00c-494c-a70b-f81c77745d21
md" #### Code 12.14"

# ╔═╡ 45fe0912-cbf9-4e82-bf7b-d60cb7aae79b
let
	pr_k = counts(trolley.response) / length(trolley.response)
	global cum_pr_k = cumsum(pr_k)
	plot(cum_pr_k, m=:o, xlab="response", ylab="cumulative proportion")
end

# ╔═╡ 1dbce450-8de8-41a6-a4ef-fd59f3663664
md" #### Code 12.15"

# ╔═╡ cd0c7f24-ec98-432c-8444-397a9bf09d91
round.(logit.(cum_pr_k), digits=2)

# ╔═╡ 67b42f0b-f991-43a9-8f9f-3a69ba540994
md" #### Code 12.16"

# ╔═╡ 068fc2dc-213a-402f-a909-0b3bcdca4dd6
@model function m12_4(R)
    # to ensure sorted cutpoints, use deltas
    Δ_cutpoints ~ filldist(Exponential(), 6)
    cutpoints = -2 .+ cumsum(Δ_cutpoints)
    R .~ OrderedLogistic(0, cutpoints)
end

# ╔═╡ 75aef58b-c087-4267-9c13-7bc1d90d1270
md" #### Code 12.18"

# ╔═╡ 13ceb9f6-46ac-473d-abe0-9b543572e5f3
begin
	Random.seed!(1)
	m12_4_ch = sample(m12_4(trolley.response), NUTS(), 1000)
	m12_4_df = DataFrame(m12_4_ch)
	PRECIS(m12_4_df)
end

# ╔═╡ 69db6156-6e68-4f4c-8ede-18968858d455
begin
	cutpoints1 = -2 .+ cumsum(mean.(eachcol(m12_4_df[!,r"Δ.*"])))
	round.(logistic.(cutpoints1), digits=3)
end

# ╔═╡ 4136174a-e7c8-4d6e-a3c8-7ed5013228fb
md" #### Code 12.19"

# ╔═╡ ba249ca4-a211-42e1-b66a-f2dd81c9b14f
md" #### Code 12.20"

# ╔═╡ dba0e917-2bd6-4fa0-b6d4-1dad4968195a
begin
	pk1 = pdf.(OrderedLogistic(0, cutpoints1), 1:7)
	round.(pk1, digits=2)
end

# ╔═╡ d0181605-1105-4fb8-8470-db8cdf3d6853
md" #### Code 12.21"

# ╔═╡ 3229c5f9-d341-4f5b-a24f-295e43437d5c
sum(pk1.*(1:7))

# ╔═╡ 696d62d7-0c39-454e-a01d-a66a54748910
md" #### Code 12.22"

# ╔═╡ 2403ce2f-7f96-4c68-9d31-80079e3daa36
begin
	pk = pdf.(OrderedLogistic(0, cutpoints1 .- 0.5), 1:7)
	round.(pk, digits=2)
end

# ╔═╡ 2a6d6a27-68b9-4722-8d51-be470942f446
md" #### Code 12.23"

# ╔═╡ 23ab19cd-6fb4-4de1-b672-daa03a1228d4
sum(pk.*(1:7))

# ╔═╡ 476d4fbd-3cc5-4e78-91f8-e8022d38ffa0
md" #### Code 12.24"

# ╔═╡ d35f0ee3-e0b7-417f-8242-6e1a5d4114d6
@model function m12_5(R, A, I, C)
    # to ensure sorted cutpoints, use deltas
    Δ_cutpoints ~ filldist(Exponential(), 6)
    cutpoints = -3 .+ cumsum(Δ_cutpoints)
    
    bA ~ Normal(0, 0.5)
    bI ~ Normal(0, 0.5)
    bC ~ Normal(0, 0.5)
    bIA ~ Normal(0, 0.5)
    bIC ~ Normal(0, 0.5)
    BI = @. bI + bIA*A + bIC*C
    phi = @. bA*A + bC*C + BI*I
    
    for i in eachindex(R)
        R[i] ~ OrderedLogistic(phi[i], cutpoints)
    end
end

# ╔═╡ 3edf24f0-be2f-461d-b121-1392618a9c67
begin
	Random.seed!(2)
	m12_5_ch = sample(m12_5(trolley.response, trolley.action, trolley.intention, trolley.contact), NUTS(), 1000)
	m12_5_df = DataFrame(m12_5_ch)
	PRECIS(m12_5_df)
end

# ╔═╡ ed58c640-6f2f-4a56-9803-8bc646a21d8f
let
	cutpoints = -2 .+ cumsum(mean.(eachcol(m12_5_df[!,r"Δ.*"])))
	round.(cutpoints; digits=3)
end

# ╔═╡ 45a8ddb4-bc79-4cfe-a3cb-07f57f441399
md" #### Code 12.25"

# ╔═╡ 8449ecae-5fd0-4c30-ba7c-582157e8d9ff
coeftab_plot(m12_5_df, pars=[:bIC, :bIA, :bC, :bI, :bA])

# ╔═╡ b268b198-44af-44da-a895-e2c245117bb5
md" #### Code 12.26"

# ╔═╡ c2a2c394-474a-4916-8e7d-5b81da7cffbe
# (not needed, in fact)


md" #### Code 12.27 & 12.28"

# ╔═╡ 3c0cdbef-6976-4b8a-8602-12706070827d
let
	p = plot(xlab="intention", ylab="probability", xlim=(0, 1), ylim=(0, 1))
	kA = 0      # value for action
	kC = 0      # value for contact
	kI = 0:1    # values for intention to calculate over
	
	rI_to_phi = (r, i) -> begin
	    BI = r.bI + r.bIA*kA + r.bIC*kC
	    r.bA*kA + r.bC*kC + BI*i
	end
	
	phi = link(m12_5_df, rI_to_phi, kI);
	p = plot(xlab="intention", ylab="probability", xlim=(0, 1), ylim=(0, 1), title="action=$kA, contact=$kC")
	for ri in 1:50
	    r = m12_5_df[ri,:]
	    cutpoints = -3 .+ cumsum(r[r"Δ.*"])
	    pk1 = cumsum(pdf.(OrderedLogistic(phi[1][ri], cutpoints), 1:6))
	    pk2 = cumsum(pdf.(OrderedLogistic(phi[2][ri], cutpoints), 1:6))
	    for i ∈ 1:6
	        plot!(kI, [pk1[i], pk2[i]], c=:gray, alpha=0.2)
	    end
	end
	p
end

# ╔═╡ 173c19d7-f51b-4260-bf8e-922e66a02721
md" #### Code 12.29"

# ╔═╡ bf96a44b-6441-4045-ae92-6613fbbbc3bb
let
	kA = 0
	kC = 1
	kI = 0:1
	
	rI_to_dist = (r, i) -> begin
	    BI = r.bI + r.bIA*kA + r.bIC*kC
	    phi = r.bA*kA + r.bC*kC + BI*i
	    cutpoints = -3 .+ cumsum(r[r"Δ.*"])
	    OrderedLogistic(phi, cutpoints)
	end
	
	Random.seed!(1)
	s = simulate(m12_5_df, rI_to_dist, kI)
	histogram(map(first, s), bar_width=0.5, label="I=0")
	histogram!(map(last, s), bar_width=0.2, label="I=1")
end

# ╔═╡ c931b144-8b63-472a-8d4c-82be650f4593
md" ## 12.4 Ordered categorical predictors"

# ╔═╡ 49bdc70e-a2d4-4238-98ac-d29b37577d80
md" #### Code 12.30"

# ╔═╡ b7ab3333-355f-42f9-b3f6-b859b3e6922d
#d = DataFrame(CSV.File("data/Trolley.csv"))
levels(trolley.edu)

# ╔═╡ 188ce9c3-9199-4768-9ad2-ed9201e563f2
md" #### Code 12.31"

# ╔═╡ 52233289-914f-4ca6-a943-9489c0158f70
edu_l = Dict{String, Int}(
    "Bachelor's Degree" => 6,
    "Elementary School" => 1,
    "Graduate Degree" => 8,
    "High School Graduate" => 4,
    "Master's Degree" => 7,
    "Middle School" => 2,
    "Some College" => 5,
    "Some High School" => 3,
)

# ╔═╡ f13534c2-2dd8-40b2-b669-7a24b5c8bb91
trolley.edu_new = map(s -> edu_l[s], trolley.edu);

# ╔═╡ eb1f70a1-42ec-4db0-b826-56bc5d4ce73b
md" #### Code 12.32"

# ╔═╡ 69afe200-63be-4a25-8b51-3a5ad992642a
begin
	Random.seed!(1805)
	delta = rand(Dirichlet(7, 2), 10)
	
	# Code 12.33
	
	h = 3
	p = plot(xlab="index", ylab="probability")
	for (idx, col) ∈ enumerate(eachcol(delta))
	    plot!(col, c=:black, alpha=0.7, lw=idx == h ? 3 : 1, m= idx == h ? :v : :o)
	end
	p
end

# ╔═╡ 6f3b01f6-b8bf-4a0e-9cc2-87c93071c9a3
md" #### Code 12.34"

# ╔═╡ 9535972c-0905-4c3a-b5b6-5656677f5d9c
# Could take 20-30 minutes...

@model function m12_6(R, action, intention, contact, E)
    delta ~ Dirichlet(7, 2)
    pushfirst!(delta, 0.0)
    
    bA ~ Normal()
    bI ~ Normal()
    bC ~ Normal()
    bE ~ Normal()

    # sum all education's deltas
    sE = sum.(map(i -> delta[1:i], E))
    phi = @. bE*sE + bA*action + bI*intention + bC*contact

    # use same cutpoints as before
    Δ_cutpoints ~ filldist(Exponential(), 6)
    cutpoints = -3 .+ cumsum(Δ_cutpoints)
    
    for i ∈ eachindex(R)
        R[i] ~ OrderedLogistic(phi[i], cutpoints)
    end
end

# ╔═╡ 7b1d63fc-bfbf-4d3e-a069-c07c7a558ee2
md" #### Code 12.35"

# ╔═╡ 859dee04-009b-4758-94ff-024a9f599509
begin
	m12_6t = m12_6(trolley.response, trolley.action, trolley.intention, trolley.contact, trolley.edu_new)
	m12_6_ch = sample(m12_6t, NUTS(), 1000);
	m12_6_df = DataFrame(m12_6_ch)
	PRECIS(m12_6_df)
end

# ╔═╡ 05c99479-f4c1-4f5f-9d92-da6b487ed00c
md" #### Code 12.36"

# ╔═╡ a1a5b1fe-1b09-4dad-8ddd-c09bbb36699f
let
	delta_labels = ["Elem","MidSch","SHS","HSG","SCol","Bach","Mast","Grad"]
	df = m12_6_df[!,r"delta.*"]
	
	corrplot(Matrix(df); seriestype=:scatter, bins=30, grid=false, ms=0.1, ma=0.8, 
	    size=(1000,800), label=delta_labels)
end

# ╔═╡ 53a9552c-f8b8-4761-8258-0062cc2ac68c
md" #### Code 12.37"

# ╔═╡ 1ac5a9e2-7e61-4f40-a8e2-c597a44bf090
trolley.edu_norm = standardize(UnitRangeTransform, Float64.(trolley.edu_new))

# ╔═╡ 6f52d2f2-008e-4991-82aa-c211cae792cd
@model function m12_7(R, action, intention, contact, edu_norm)
    bA ~ Normal()
    bI ~ Normal()
    bC ~ Normal()
    bE ~ Normal()

    phi = @. bE*edu_norm + bA*action + bI*intention + bC*contact

    # use same cutpoints as before
    Δ_cutpoints ~ filldist(Exponential(), 6)
    cutpoints = -3 .+ cumsum(Δ_cutpoints)
    
    for i ∈ eachindex(R)
        R[i] ~ OrderedLogistic(phi[i], cutpoints)
    end
end

# ╔═╡ dffc67be-7ab2-49ba-a584-7877faed20d6
begin
	model = m12_7(trolley.response, trolley.action, trolley.intention, trolley.contact, trolley.edu_norm)
	m12_7_ch = sample(model, NUTS(), 1000)
	m12_7_df = DataFrame(m12_7_ch)
	PRECIS(m12_7_df)
end

# ╔═╡ ce0ee122-a7ed-417d-abaa-ae25aeeaabb4
import BASE: rand

# ╔═╡ 9f5a39bc-9e1c-4f08-ba50-0cf87122358b
rand(d::ZIPoisson, N::Int) = map(_->rand(d), 1:N)

# ╔═╡ 0ebe4814-ce59-4354-87de-e80e540e355c
function rand(d::ZIPoisson)
    rand() <= d.w ? 0 : rand(Poisson(d.λ))
end

# ╔═╡ Cell order:
# ╠═147708bd-6ded-4316-aca0-1f6524ba3c53
# ╠═2b447c25-cc6e-4b40-a18a-fcab8eb64d8e
# ╠═c9ea9876-b966-41b4-ac02-a6b1040e7119
# ╠═2c26b000-36ee-473a-a7c5-d6a746096ec9
# ╟─7c87ccc8-8b76-45f1-b10f-b04bba307ce9
# ╟─5be273e4-faae-4f06-934c-a1da019e8b89
# ╠═1651ab99-e88a-40b7-9092-daf7f31e5604
# ╠═fdd7dad1-facd-4ce3-9a72-20ce2c4cbe6a
# ╟─86d02852-70e2-4fcc-bc8c-8d2b881fa248
# ╠═f2af75cc-361e-4a45-a1b4-14a4ce404fa3
# ╠═8d1b576f-e2a4-423d-80cc-0e69f6a376d4
# ╠═9ce82e0e-bdcf-426c-b87e-32346c658b87
# ╟─b339874c-e364-4f8c-bd62-5d66b175f0be
# ╠═50661b7a-3484-4c52-b87c-4f44f119e164
# ╟─4335575b-2ce8-41f0-a7ed-295789a5b723
# ╠═c8959699-8886-423a-b994-3180c639de1b
# ╟─7959fa22-560e-44d0-96ee-d2e7bf51cbf9
# ╠═ce0ee122-a7ed-417d-abaa-ae25aeeaabb4
# ╠═c27a8622-edca-4f19-91d0-1682ebc9cc7d
# ╠═c799de28-efa8-4466-812f-ac3c79b5d835
# ╟─7a9ad8dc-019d-4a32-b8b3-ee5a2703e609
# ╠═5e3a0ec2-5066-4ab6-881d-af6113d23974
# ╠═29167add-b592-4a1e-8742-88fa3264286d
# ╠═cc4b5a4e-c2c6-4295-b4ff-8f73e84fcd88
# ╟─512a3b1f-2e79-4777-b976-6d310725e9f8
# ╟─91de0472-61b8-47f4-babb-e07539fa3a55
# ╠═15d28187-c155-4332-bde4-90c2d40b7a36
# ╟─dc4292b2-5aa8-4413-b4a3-91456f1d1d2a
# ╠═3288b258-640b-414f-90d4-0bc1517aa506
# ╠═45e9dc3c-bd28-4eaf-a56a-d12f4804f341
# ╠═ab414dc7-7f31-4b0a-9d65-1268940fbc50
# ╠═0ebe4814-ce59-4354-87de-e80e540e355c
# ╠═9f5a39bc-9e1c-4f08-ba50-0cf87122358b
# ╠═73d535fd-227a-4c4b-a0ca-ed926ecf86c6
# ╠═cb63530d-4041-4a53-92c7-73a5c5705564
# ╟─3b4f512c-7fc3-4eda-90a3-57d7ed313f87
# ╠═719b6fe8-03e0-4e1d-8be8-6d642f78e86d
# ╟─e3ea5e44-b189-4c2e-b7d0-5bdc26477ef8
# ╟─8c89b94a-3250-4abf-8eda-35fd2b42f708
# ╠═a6bcc9c6-4be7-4390-9fa3-6a124ed8afc1
# ╟─cbf2db2b-5e72-4ad1-9fc3-0944c17db2df
# ╠═c1943565-d495-480d-86ce-ea6f4e760b53
# ╟─ee421d11-e00c-494c-a70b-f81c77745d21
# ╠═45fe0912-cbf9-4e82-bf7b-d60cb7aae79b
# ╠═1dbce450-8de8-41a6-a4ef-fd59f3663664
# ╠═cd0c7f24-ec98-432c-8444-397a9bf09d91
# ╟─67b42f0b-f991-43a9-8f9f-3a69ba540994
# ╠═068fc2dc-213a-402f-a909-0b3bcdca4dd6
# ╟─75aef58b-c087-4267-9c13-7bc1d90d1270
# ╠═13ceb9f6-46ac-473d-abe0-9b543572e5f3
# ╠═69db6156-6e68-4f4c-8ede-18968858d455
# ╟─4136174a-e7c8-4d6e-a3c8-7ed5013228fb
# ╟─ba249ca4-a211-42e1-b66a-f2dd81c9b14f
# ╠═dba0e917-2bd6-4fa0-b6d4-1dad4968195a
# ╟─d0181605-1105-4fb8-8470-db8cdf3d6853
# ╠═3229c5f9-d341-4f5b-a24f-295e43437d5c
# ╟─696d62d7-0c39-454e-a01d-a66a54748910
# ╠═2403ce2f-7f96-4c68-9d31-80079e3daa36
# ╟─2a6d6a27-68b9-4722-8d51-be470942f446
# ╠═23ab19cd-6fb4-4de1-b672-daa03a1228d4
# ╟─476d4fbd-3cc5-4e78-91f8-e8022d38ffa0
# ╠═d35f0ee3-e0b7-417f-8242-6e1a5d4114d6
# ╠═3edf24f0-be2f-461d-b121-1392618a9c67
# ╠═ed58c640-6f2f-4a56-9803-8bc646a21d8f
# ╠═45a8ddb4-bc79-4cfe-a3cb-07f57f441399
# ╠═8449ecae-5fd0-4c30-ba7c-582157e8d9ff
# ╠═b268b198-44af-44da-a895-e2c245117bb5
# ╟─c2a2c394-474a-4916-8e7d-5b81da7cffbe
# ╠═3c0cdbef-6976-4b8a-8602-12706070827d
# ╟─173c19d7-f51b-4260-bf8e-922e66a02721
# ╠═bf96a44b-6441-4045-ae92-6613fbbbc3bb
# ╟─c931b144-8b63-472a-8d4c-82be650f4593
# ╟─49bdc70e-a2d4-4238-98ac-d29b37577d80
# ╠═b7ab3333-355f-42f9-b3f6-b859b3e6922d
# ╟─188ce9c3-9199-4768-9ad2-ed9201e563f2
# ╠═52233289-914f-4ca6-a943-9489c0158f70
# ╠═f13534c2-2dd8-40b2-b669-7a24b5c8bb91
# ╟─eb1f70a1-42ec-4db0-b826-56bc5d4ce73b
# ╠═69afe200-63be-4a25-8b51-3a5ad992642a
# ╟─6f3b01f6-b8bf-4a0e-9cc2-87c93071c9a3
# ╠═9535972c-0905-4c3a-b5b6-5656677f5d9c
# ╟─7b1d63fc-bfbf-4d3e-a069-c07c7a558ee2
# ╠═859dee04-009b-4758-94ff-024a9f599509
# ╟─05c99479-f4c1-4f5f-9d92-da6b487ed00c
# ╠═a1a5b1fe-1b09-4dad-8ddd-c09bbb36699f
# ╟─53a9552c-f8b8-4761-8258-0062cc2ac68c
# ╠═1ac5a9e2-7e61-4f40-a8e2-c597a44bf090
# ╠═6f52d2f2-008e-4991-82aa-c211cae792cd
# ╠═dffc67be-7ab2-49ba-a584-7877faed20d6
