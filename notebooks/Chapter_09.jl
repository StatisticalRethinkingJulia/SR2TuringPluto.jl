### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 4b348c51-c844-4e64-85e4-a4f7ea952fc9
using Pkg, DrWatson

# ╔═╡ 75ce54b9-b801-412a-b70e-ab5222df137f
begin
	using Random
	using StatsBase
	using Distributions
	using StatsPlots
	using StatsFuns
	using Logging
	
	using CSV
	using DataFrames
	using Optim
	
	using MCMCChains
	using Optim
	using Turing
	using StatisticalRethinking
end

# ╔═╡ e803d172-2b8a-43b4-b99d-307e64cf29d5
md"# Chapter 9. Markov Chain Monte Carlo."

# ╔═╡ 77a70403-fa95-4b13-b886-415d1a91f15a
begin
	default(labels=false)
	Logging.disable_logging(Logging.Warn);
end

# ╔═╡ f9ff23ea-fdf2-4d48-bb84-98f40db8603e
md"## 9.1 Good King Markov and his island kingdom."

# ╔═╡ 821bcd65-4ab5-43e8-837e-8bdf96f24edd
md"### Code 9.1"

# ╔═╡ 0ccba2cd-c10c-4713-80e1-77e4a3e9c85c
begin
	Random.seed!(1)
	num_weeks = 10^5
	positions = []
	current = 10
end;

# ╔═╡ b71b59cb-b9ed-4a18-a713-400d8fb48517
for i ∈ 1:num_weeks
    # record current position
    push!(positions, current)
    # flip coin to generate proposal
    proposal = current + sample([-1, 1])
    # handle loops around
    proposal < 1 && (proposal = 10)
    proposal > 10 && (proposal = 1)
    # move?
    prob_move = proposal / current
    rand() < prob_move && (current = proposal)
end

# ╔═╡ 62852231-aebe-43bc-a591-ff3c79651739
md"### Code 9.2"

# ╔═╡ 2f603bc2-938d-45b4-a9ed-60daa5847af5
scatter(positions[1:100], xlab="week", ylab="island")

# ╔═╡ 3f745e69-690e-4ab9-934f-b05f24f20a98
md"### Code 9.3"

# ╔═╡ a88970a9-341c-4e92-891b-6057dbac1e01
histogram(positions, xlab="island", ylab="number of weeks")

# ╔═╡ e7569635-07ab-4b69-8014-f26d42955fce
md"## 9.2 Metropolis algorithms"

# ╔═╡ d4991ad9-4341-47e8-aee4-9c35cbe28e35
md"### Code 9.4"

# ╔═╡ 9322d404-9ebd-4dd8-985d-8322f4f77571
begin
	D = 10
	T = 1000
	Y = rand(MvNormal(zeros(D), ones(D)), T)
	Rd = sqrt.(sum.(eachcol(Y.^2)))
	density(Rd)
end

# ╔═╡ 7d94a7f1-7ce5-4f6a-ba72-0a8a9cde44e9
md"## 9.3 Hamiltonian Monte Carlo"

# ╔═╡ 19d13362-93a0-4a93-b3bb-6e30d3d934a3
md"### Code 9.5"

# ╔═╡ 29a1a866-522b-4354-866e-9f121ba378b9
begin
	Random.seed!(7)
	
	x = rand(Normal(), 50)
	y = rand(Normal(), 50)
	x = standardize(ZScoreTransform, x)
	y = standardize(ZScoreTransform, y);
end;

# ╔═╡ e3cc237b-9997-41b6-8298-1b79e4444662
function U(q::Vector{Float64}; a=0, b=1, k=0, d=1)::Float64
    μy, μx = q
    U = sum(normlogpdf.(μy, 1, y)) + sum(normlogpdf.(μx, 1, x)) 
    U += normlogpdf(a, b, μy) + normlogpdf(k, d, μx)
    -U
end

# ╔═╡ cb41155d-43a6-4fef-932f-13cabcdf907a
md"### Code 9.6"

# ╔═╡ 7d29cbaf-2b27-4b4e-8dee-74a049408fba
function ∇U(q::Vector{Float64}; a=0, b=1, k=0, d=1)::Vector{Float64}
    μy, μx = q
    G₁ = sum(y .- μy) + (a - μy) / b^2  # ∂U/∂μy
    G₂ = sum(x .- μx) + (k - μx) / d^2  # ∂U/∂μx
    [-G₁, -G₂]
end

# ╔═╡ 9d674865-3fed-48a6-a970-d997ae9d7443
md"### Codes 9.8 - 9.10 (before 9.7 to define HMC2 function)"

# ╔═╡ aa9be028-c3dc-4e4a-9405-86be642540ca
function HMC2(U, ∇U, ϵ::Float64, L::Int, current_q::Vector{Float64})
    q = current_q
    p = rand(Normal(), length(q))  # random flick - p is momentum
    current_p = p
    
    # make a half step for momentum at the beginning
    p -= ϵ .* ∇U(q) ./ 2
    
    # initialize bookkeeping - saves trajectory
    qtraj = [q]
    ptraj = [p]
    
    # Alternate full steps for position and momentum
    for i ∈ 1:L
        q += @. ϵ * p  # full step for the position
        # make a full step for the momentum except at the end of trajectory
        if i != L
            p -= ϵ * ∇U(q)
            push!(ptraj, p)
        end
        push!(qtraj, q)
    end
    
    # Make a half step for momentum at the end
    p -= ϵ * ∇U(q) / 2
    push!(ptraj, p)
    
    # negate momentum at the end of trajectory to make the proposal symmetric
    p = -p
    
    # evaluate potential and kinetic energies at the start and the end of trajectory
    current_U = U(current_q)
    current_K = sum(current_p.^2)/2
    proposed_U = U(q)
    proposed_K = sum(p.^2)/2
    
    # accept or reject the state at the end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    accept = (rand() < exp(current_U - proposed_U + current_K - proposed_K))

    if accept
        current_q = q
    end
    
    (q=current_q, traj=qtraj, ptraj=ptraj, accept=accept)
end

# ╔═╡ e61d81ed-761a-450b-9ea0-1398a9ac8e22
md"### Code 9.7"

# ╔═╡ a74b283d-febc-40b1-9654-3d63b7d014b2
begin
	Random.seed!(1)
	Q = (q=[-0.1, 0.2],)
	pr = 0.3
	step1 = 0.03
	L = 11
	n_samples = 4
	p = scatter([Q.q[1]], [Q.q[2]], xlab="μx", ylab="μy")


	for i ∈ 1:n_samples
	    Q = HMC2(U, ∇U, step1, L, Q.q)
	    if n_samples < 10 
	        cx, cy = [], []
	        for j ∈ 1:L
	            K0 = sum(Q.ptraj[j].^2)/2
	            plot!(
	                [Q.traj[j][1], Q.traj[j+1][1]],
	                [Q.traj[j][2], Q.traj[j+1][2]],
	                lw=1+2*K0,
	                c=:black,
	                alpha=0.5
	            )
	            push!(cx, Q.traj[j+1][1])
	            push!(cy, Q.traj[j+1][2])
	        end
	        scatter!(cx, cy, c=:white, ms=3)
	    end
	    scatter!([Q.q[1]], [Q.q[2]], shape=(Q.accept ? :circle : :rect), c=:blue)
	end
	p
end

# ╔═╡ 4dc8aa1b-7544-4a80-9109-a2c9040d9321
md"## 9.4 Easy HMC: ulam"

# ╔═╡ e6a3d465-57d2-40aa-b71d-8cca5db234d9
md"### Code 9.11"

# ╔═╡ ff957f53-dbc7-4885-9797-61d10b8dac15
begin
	d = CSV.read(sr_datadir("rugged.csv"), DataFrame)
	dd = d[completecases(d, :rgdppc_2000),:]
	dd[:,:log_gdp] = log.(dd.rgdppc_2000);
	dd[:,:log_gdp_std] = dd.log_gdp / mean(dd.log_gdp)
	dd[:,:rugged_std] = dd.rugged / maximum(dd.rugged)
	dd[:,:cid] = @. ifelse(dd.cont_africa == 1, 1, 2);
end;

# ╔═╡ 06d5c8b7-7a87-4ddb-9c30-d6d623a71814
md"### Code 9.12"

# ╔═╡ e29edd1a-890d-4e21-a221-d4311bd845d3
r̄ = mean(dd.rugged_std);

# ╔═╡ 200f0e88-b30f-4aac-a733-a25f5a52c5d8
@model function model_m8_3(rugged_std, cid,  log_gdp_std)
    σ ~ Exponential()
    a ~ MvNormal([1, 1], 0.1)
    b ~ MvNormal([0, 0], 0.3)
    μ = @. a[cid] + b[cid] * (rugged_std - r̄)
    log_gdp_std ~ MvNormal(μ, σ)
end

# ╔═╡ ceb3dbdd-a7d3-4eac-970b-3b10c9869974
m8_3 = optimize(model_m8_3(dd.rugged_std, dd.cid,  dd.log_gdp_std), MAP())

# ╔═╡ f67d530a-e56d-408c-8f21-d4e37e8dab42
md"### Code 9.13"

# ╔═╡ 60a185df-0656-4a46-b301-388407073ab6
md"#### For Turing this is not needed"

# ╔═╡ 9d6afd05-6570-44e0-8b49-5e4bc7e8a76f
begin
	dat_slim = dd[!,[:log_gdp_std, :rugged_std, :cid]]
	describe(dat_slim)
end

# ╔═╡ cfe67e9f-41bb-419f-96ba-282926de31a2
md"### Code 9.14"

# ╔═╡ 8409b9f4-7fa5-4f38-9b38-ae34e4c60e83
@model function model_m9_1(rugged_std, cid,  log_gdp_std)
    σ ~ Exponential()
    a ~ MvNormal([1, 1], 0.1)
    b ~ MvNormal([0, 0], 0.3)
    μ = @. a[cid] + b[cid] * (rugged_std - r̄)
    log_gdp_std ~ MvNormal(μ, σ)
end

# ╔═╡ 4e94cf59-81eb-4b99-b048-502f54faee15
md"#### One chain will be produced by default"

# ╔═╡ 76782d05-9e32-424c-a206-8081e0fdd146
m9_1 = sample(model_m8_3(dd.rugged_std, dd.cid,  dd.log_gdp_std), NUTS(), 1000);

# ╔═╡ e8dd780e-b94d-413c-9b0b-0476817b2094
md"### Code 9.15"

# ╔═╡ b80e5ff3-445a-465c-9665-697f5b0b897e
describe(DataFrame(m9_1))

# ╔═╡ b6057787-3948-482c-be7e-69c4231386b6
md"### Code 9.16"

# ╔═╡ d438bddd-9138-410e-ab98-40f7c90e780f
md"#### For this to use multiple cores, julia has to be started with `--threads 4` parameter, otherwise chains will be sampled sequentially"

# ╔═╡ 3f4a3b7c-2db7-4d7b-800a-c14968a02855
m9_1_4 = sample(model_m8_3(dd.rugged_std, dd.cid,  dd.log_gdp_std), NUTS(),
	, 1000);

# ╔═╡ 2d2aa3ef-8fb7-4c69-8869-42e1939dc950
md"### Code 9.17"

# ╔═╡ 61b5e6a8-3baa-4101-9a21-5d49af43bfe7
md"#### This shows combined chains statistics. To get information about individual chains, use `m9_1[:,:,1]`"

# ╔═╡ 35f25d1e-234b-481f-9303-3c5973ab1228
m9_1

# ╔═╡ 4e20ad94-3011-4b89-b29e-8a793e733e03
md"### Code 9.18"

# ╔═╡ ed026592-6cb3-40a3-9f0a-565e14038114
describe(DataFrame(m9_1[:,:,1]))

# ╔═╡ b25c0b7a-4469-4945-be28-bfb2c22103ca
md"### Code 9.19"

# ╔═╡ 8a0295cc-4f0c-475a-a8c6-ccf2b99ed00f
@df DataFrame(m9_1) corrplot(cols(1:5), seriestype=:scatter, ms=0.2, size=(950, 800), bins=30, grid=false)

# ╔═╡ 580db5a8-550a-42d1-9aef-c95f4099a3fa
md"### Code 9.20"

# ╔═╡ 804db079-331c-4f76-b62f-bbc66911b5f2
plot(m9_1_4)

# ╔═╡ 534f0c53-da8b-4da5-81d4-c0ab3ceb270b
md"### Code 9.21"

# ╔═╡ c634d1a4-5f08-4bea-b5ab-acc8ab3af668
plot(m9_1)

# ╔═╡ 4535c71f-be7c-41e8-940b-57289c83e509
md"## 9.5 Care and feeding of your Markov chain."

# ╔═╡ f39cd55f-ed8f-487e-885e-f1e16b78d9dd
md"### Codes 9.22 - 9.23"

# ╔═╡ cd971faa-315a-4cea-a11a-0cb28f3e3997
let
	# To make it diverting with Turing, it was needed to increase exp() argument.
	Random.seed!(1)
	y = [-1., 1.]

	@model function model_m9_2(y)
	    α ~ Normal(0, 1000)
	    σ ~ Exponential(1/0.0001)
	    y ~ Normal(α, σ)
	end

	global m9_2 = sample(model_m9_2(y), NUTS(), 1000)
	m9_2_df = DataFrame(m9_2)
	describe(m9_2_df)
end

# ╔═╡ a9febef3-06d0-4552-8590-2f165ce36547
md"### Code 9.23"

# ╔═╡ 7b2588d2-d86d-44fb-9895-d6452875bc47
plot(
    traceplot(m9_2),
    histogram(m9_2),
    size=(900, 500)
)

# ╔═╡ cc97f4d9-ad00-4d2c-b520-2071b9cf86f6
md"### Code 9.24"

# ╔═╡ 631bb03b-00a5-466a-be0d-36f86a9721aa
Random.seed!(2)

# ╔═╡ 95dc66d6-660a-4b0c-b89d-487a22dfec82
@model function model_m9_3(y)
    α ~ Normal(1, 10)
    σ ~ Exponential(1)
    y ~ Normal(α, σ)
end

# ╔═╡ c5570bf1-c04d-41e3-b8ae-6e6dac45cb70
begin
	m9_3 = sample(model_m9_3(y), NUTS(), 1000)
	m9_3_df = DataFrame(m9_3)
	describe(m9_3_df)
end

# ╔═╡ 8f945f6c-2561-433b-aedb-3448d2c4b419
ess_rhat(m9_3)

# ╔═╡ ecd23fa3-fdd2-4ffb-a0e9-ad1feee04d22
md"### Code 9.25 - 9.26"

# ╔═╡ f10e310f-165b-4b4b-b3e0-66af6430ba4f
let
	Random.seed!(41)
	y = rand(Normal(), 100)

	Random.seed!(384)

	@model function model_m9_4(y)
	    a1 ~ Normal(0, 1000)
	    a2 ~ Normal(0, 1000)
	    σ ~ Exponential(1)
	    μ = a1 + a2
	    y ~ Normal(μ, σ)
	end

	global m9_4 = sample(model_m9_4(y), NUTS(), 1000)
	m9_4_df = DataFrame(m9_4)
	describe(m9_4_df)

end

# ╔═╡ 4c09eda3-8f4c-4fa2-91ea-99b92cfbb422
ess_rhat(m9_4)

# ╔═╡ 1f2c5378-bbcd-4fb7-9130-9c6cebb4c6d4
plot(m9_4)

# ╔═╡ 517d31ba-5843-40a7-ad06-96e7a8ea8845
md"### Code 9.27"

# ╔═╡ 1cc0a6c1-3fcf-4df3-a8d9-ba3454c0a37b
Random.seed!(384)

# ╔═╡ 9a76593b-c9a9-45cc-aca5-9877432baffd
@model function model_m9_5(y)
    a1 ~ Normal(0, 10)
    a2 ~ Normal(0, 10)
    σ ~ Exponential(1)
    μ = a1 + a2
    y ~ Normal(μ, σ)
end

# ╔═╡ 6a7ea8be-bea2-41c4-89ab-9a6bbcfba116
begin
	m9_5 = sample(model_m9_5(y), NUTS(), 1000)
	m9_5_df = DataFrame(m9_5)
	describe(m9_5_df)
end

# ╔═╡ 942e3f19-a5c1-430e-a3be-097c53bc7394
ess_rhat(m9_5)

# ╔═╡ 6caa4b75-9ac1-496b-a4b2-38b72b203e59
plot(m9_5)

# ╔═╡ Cell order:
# ╟─e803d172-2b8a-43b4-b99d-307e64cf29d5
# ╠═4b348c51-c844-4e64-85e4-a4f7ea952fc9
# ╠═75ce54b9-b801-412a-b70e-ab5222df137f
# ╠═77a70403-fa95-4b13-b886-415d1a91f15a
# ╟─f9ff23ea-fdf2-4d48-bb84-98f40db8603e
# ╟─821bcd65-4ab5-43e8-837e-8bdf96f24edd
# ╠═0ccba2cd-c10c-4713-80e1-77e4a3e9c85c
# ╠═b71b59cb-b9ed-4a18-a713-400d8fb48517
# ╟─62852231-aebe-43bc-a591-ff3c79651739
# ╠═2f603bc2-938d-45b4-a9ed-60daa5847af5
# ╟─3f745e69-690e-4ab9-934f-b05f24f20a98
# ╠═a88970a9-341c-4e92-891b-6057dbac1e01
# ╟─e7569635-07ab-4b69-8014-f26d42955fce
# ╟─d4991ad9-4341-47e8-aee4-9c35cbe28e35
# ╠═9322d404-9ebd-4dd8-985d-8322f4f77571
# ╟─7d94a7f1-7ce5-4f6a-ba72-0a8a9cde44e9
# ╟─19d13362-93a0-4a93-b3bb-6e30d3d934a3
# ╠═29a1a866-522b-4354-866e-9f121ba378b9
# ╠═e3cc237b-9997-41b6-8298-1b79e4444662
# ╟─cb41155d-43a6-4fef-932f-13cabcdf907a
# ╠═7d29cbaf-2b27-4b4e-8dee-74a049408fba
# ╟─9d674865-3fed-48a6-a970-d997ae9d7443
# ╠═aa9be028-c3dc-4e4a-9405-86be642540ca
# ╟─e61d81ed-761a-450b-9ea0-1398a9ac8e22
# ╠═a74b283d-febc-40b1-9654-3d63b7d014b2
# ╟─4dc8aa1b-7544-4a80-9109-a2c9040d9321
# ╟─e6a3d465-57d2-40aa-b71d-8cca5db234d9
# ╠═ff957f53-dbc7-4885-9797-61d10b8dac15
# ╟─06d5c8b7-7a87-4ddb-9c30-d6d623a71814
# ╠═e29edd1a-890d-4e21-a221-d4311bd845d3
# ╠═200f0e88-b30f-4aac-a733-a25f5a52c5d8
# ╠═ceb3dbdd-a7d3-4eac-970b-3b10c9869974
# ╟─f67d530a-e56d-408c-8f21-d4e37e8dab42
# ╟─60a185df-0656-4a46-b301-388407073ab6
# ╠═9d6afd05-6570-44e0-8b49-5e4bc7e8a76f
# ╟─cfe67e9f-41bb-419f-96ba-282926de31a2
# ╠═8409b9f4-7fa5-4f38-9b38-ae34e4c60e83
# ╟─4e94cf59-81eb-4b99-b048-502f54faee15
# ╠═76782d05-9e32-424c-a206-8081e0fdd146
# ╟─e8dd780e-b94d-413c-9b0b-0476817b2094
# ╠═b80e5ff3-445a-465c-9665-697f5b0b897e
# ╟─b6057787-3948-482c-be7e-69c4231386b6
# ╟─d438bddd-9138-410e-ab98-40f7c90e780f
# ╠═3f4a3b7c-2db7-4d7b-800a-c14968a02855
# ╟─2d2aa3ef-8fb7-4c69-8869-42e1939dc950
# ╟─61b5e6a8-3baa-4101-9a21-5d49af43bfe7
# ╠═35f25d1e-234b-481f-9303-3c5973ab1228
# ╟─4e20ad94-3011-4b89-b29e-8a793e733e03
# ╠═ed026592-6cb3-40a3-9f0a-565e14038114
# ╟─b25c0b7a-4469-4945-be28-bfb2c22103ca
# ╠═8a0295cc-4f0c-475a-a8c6-ccf2b99ed00f
# ╟─580db5a8-550a-42d1-9aef-c95f4099a3fa
# ╠═804db079-331c-4f76-b62f-bbc66911b5f2
# ╟─534f0c53-da8b-4da5-81d4-c0ab3ceb270b
# ╠═c634d1a4-5f08-4bea-b5ab-acc8ab3af668
# ╟─4535c71f-be7c-41e8-940b-57289c83e509
# ╟─f39cd55f-ed8f-487e-885e-f1e16b78d9dd
# ╠═cd971faa-315a-4cea-a11a-0cb28f3e3997
# ╟─a9febef3-06d0-4552-8590-2f165ce36547
# ╠═7b2588d2-d86d-44fb-9895-d6452875bc47
# ╟─cc97f4d9-ad00-4d2c-b520-2071b9cf86f6
# ╠═631bb03b-00a5-466a-be0d-36f86a9721aa
# ╠═95dc66d6-660a-4b0c-b89d-487a22dfec82
# ╠═c5570bf1-c04d-41e3-b8ae-6e6dac45cb70
# ╠═8f945f6c-2561-433b-aedb-3448d2c4b419
# ╟─ecd23fa3-fdd2-4ffb-a0e9-ad1feee04d22
# ╠═f10e310f-165b-4b4b-b3e0-66af6430ba4f
# ╠═4c09eda3-8f4c-4fa2-91ea-99b92cfbb422
# ╠═1f2c5378-bbcd-4fb7-9130-9c6cebb4c6d4
# ╟─517d31ba-5843-40a7-ad06-96e7a8ea8845
# ╠═1cc0a6c1-3fcf-4df3-a8d9-ba3454c0a37b
# ╠═9a76593b-c9a9-45cc-aca5-9877432baffd
# ╠═6a7ea8be-bea2-41c4-89ab-9a6bbcfba116
# ╠═942e3f19-a5c1-430e-a3be-097c53bc7394
# ╠═6caa4b75-9ac1-496b-a4b2-38b72b203e59
