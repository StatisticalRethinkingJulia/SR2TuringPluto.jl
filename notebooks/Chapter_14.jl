### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 4458903f-4717-4603-a707-4d2fd587470b
using Pkg, DrWatson

# ╔═╡ d577b8c8-3670-497a-a6fe-41208a8e8501
begin
	using Turing
	using Turing
	using DataFrames
	using CSV
	using Random
	using Dagitty
	using Distributions
	using StatisticalRethinking
	using StatisticalRethinking: link
	using StatisticalRethinkingPlots
	using StatsPlots
	using StatsBase
	using Logging
	using LinearAlgebra
end

# ╔═╡ 95965ab2-f250-479f-9d5d-e478dba1fe35
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

# ╔═╡ 90028541-7308-476f-b51b-850cd6ad39c4
begin
	default(label=false);
	Logging.disable_logging(Logging.Warn)
end;

# ╔═╡ 1217e4ab-8e97-4b73-9493-a58c7d8165fe
md" ## 14.1 Varying slopes by construction."

# ╔═╡ 100c3470-b7e4-4890-9347-851b78f5202f
md" #### Code 14.1 - 14.4"

# ╔═╡ 464ddfd5-d06b-4615-bd23-e63e03a29f1e
md"

!!! note

Julia has similar `column-first` matix order."

# ╔═╡ fe81bd6a-c7bf-437d-8ca6-f21f9b45b41f
reshape(1:4, (2,2))

# ╔═╡ b6c8f7af-8e39-4792-a62c-88a41876345d
md" #### Code 14.5 - 14.8"

# ╔═╡ 8ee2eec3-3214-48d1-8604-a54ba1a67baf
let
	a = 3.5    # average morning wait time
	b = -1     # average difference afternoon wait time
	σ_a = 1    # std dev in intercepts
	σ_b = 0.5  # std dev in slopes
	ρ = -0.7;  # correlation between intercepts and slopes
	
	global μ = [a, b]
	cov_ab = σ_a * σ_b * ρ
	global Σ₂ = [[σ_a^2, cov_ab] [cov_ab, σ_b^2]]
	sigmas = [σ_a, σ_b]
	Ρ = [[1, ρ] [ρ, 1]]
	
	global Σ₃ = Diagonal(sigmas) * Ρ * Diagonal(sigmas)
	N_cafes = 20
	Random.seed!(5)
	vary_effect = rand(MvNormal(μ, Σ₃), N_cafes)
	global a_cafe = vary_effect[1,:]
	global b_cafe = vary_effect[2,:];
end

# ╔═╡ 573253a5-10d4-458e-9451-89c763658810
md" #### Code 14.9"

# ╔═╡ fa9dff6f-3182-4093-bbfd-8deea404eab0
let
	p = scatter(a_cafe, b_cafe, xlab="intercepts (a_cafe)", ylab="slopes (b_cafe)")
	
	d = acos(Σ₂[1,2])
	chi = Chisq(2)
	
	for l ∈ (0.1, 0.3, 0.5, 0.8, 0.99)
	    scale = sqrt(quantile(chi, l))
	    xₜ(t) = scale*Σ₂[1,1]*cos(t + d/2) + μ[1]
	    yₜ(t) = scale*Σ₂[2,2]*cos(t - d/2) + μ[2]
	
	    plot!(xₜ, yₜ, 0, 2π, c=:black, alpha=0.3)
	end
	p
end

# ╔═╡ 2367f3cf-c568-4f6c-82c6-b7ae42408404
md" #### Code 14.10"

# ╔═╡ 8185285a-1640-45e9-8c80-b17b6c960804
let
	Random.seed!(1)
	N_cafes = 20
	N_visits = 10
	
	afternoon = repeat(0:1, N_visits*N_cafes ÷ 2)
	cafe_id = repeat(1:N_cafes, inner=N_visits)
	μ = a_cafe[cafe_id] + b_cafe[cafe_id] .* afternoon
	σ = 0.5
	wait = rand.(Normal.(μ, σ))
	global cafes = DataFrame(cafe=cafe_id, afternoon=afternoon, wait=wait);
end

# ╔═╡ 40634f53-8104-4791-877c-2d498cab9570
md" #### Code 14.11"

# ╔═╡ 885a0261-3541-4389-beb7-231c77e2dfcc
let
	R = rand(LKJ(2, 2), 10^4);
	density(getindex.(R, 2), xlab="Correlation")
end

# ╔═╡ 6837fff8-0c29-4ba1-9b3b-b21c23d861eb
md" #### Code 14.12"

# ╔═╡ 1f30d601-4e00-44dd-a6ab-363c431ca875
@model function m14_1(cafe, afternoon, wait)
    a ~ Normal(5, 2)
    b ~ Normal(-1, 0.5)
    σ_cafe ~ filldist(Exponential(), 2)
    Rho ~ LKJ(2, 2)
    # build sigma matrix manually, to avoid numerical errors
#     (σ₁, σ₂) = σ_cafe
#     sc = [[σ₁^2, σ₁*σ₂] [σ₁*σ₂, σ₂^2]]
#     Σ = Rho .* sc
    # the same as above, but shorter and generic
    Σ = (σ_cafe .* σ_cafe') .* Rho
    ab ~ filldist(MvNormal([a,b], Σ), 20)
    a = ab[1,cafe]
    b = ab[2,cafe]
    μ = @. a + b * afternoon
    σ ~ Exponential()
    for i ∈ eachindex(wait)
        wait[i] ~ Normal(μ[i], σ)
    end
end

# ╔═╡ de06a783-c2ff-4f45-8a1e-02977ad31e26
begin
	Random.seed!(1)
	m14_1_ch = sample(m14_1(cafes.cafe, cafes.afternoon, cafes.wait), NUTS(), 1000)
	m14_1_df = DataFrame(m14_1_ch);
end

# ╔═╡ 30be727e-86e6-41cf-be58-e0b1ffb372c2
md" #### Code 14.13"

# ╔═╡ 5ca07e4a-bf45-4e3c-9a5b-3b9b7af38020
let
	density(m14_1_df."Rho[1,2]", lab="posterior", lw=2)
	R = rand(LKJ(2, 2), 10^4);
	density!(getindex.(R, 2), lab="prior", ls=:dash, lw=2)
	plot!(xlab="correlation", ylab="Density")
end

# ╔═╡ f95ec79a-0263-46e4-b82e-a274cb7ee7fb
md" #### Code 14.14"

# ╔═╡ 81f57cdb-36f8-43bb-86a5-19bac5cd88b1
md"

!!! note

Plot differs from presented in the book due to different seed in data generation."

# ╔═╡ 78fe51f6-beea-4155-905d-2d5e6e7d66d7
let
	N_cafes = 20
	gb = groupby(cafes[cafes.afternoon .== 0,:], :cafe)
	global a1 = combine(gb, :wait => mean).wait_mean
	
	gb = groupby(cafes[cafes.afternoon .== 1,:], :cafe)
	global b1 = combine(gb, :wait => mean).wait_mean .- a1
	
	global a2 = [mean(m14_1_df[:, "ab[1,$i]"]) for i ∈ 1:N_cafes]
	global b2 = [mean(m14_1_df[:, "ab[2,$i]"]) for i ∈ 1:N_cafes]
	
	xlim = extrema(a1) .+ (-0.1, 0.1)
	ylim = extrema(b1) .+ (-0.1, 0.1)
	
	global p_cafe = scatter(a1, b1, xlab="intercept", ylab="slope", xlim=xlim, ylim=ylim)
	
	scatter!(a2, b2, mc=:white)
	
	for (x1,y1,x2,y2) ∈ zip(a1, b1, a2, b2)
	    plot!([x1,x2], [y1,y2], c=:black)
	end
	p_cafe
end

# ╔═╡ 4ea3a429-e457-48a1-a66a-ea66ff21d22c
md" #### Code 14.15"

# ╔═╡ 8bb5afba-3f59-4776-9a7b-8df9fed58a90
let
	# posterior mean
	ρ = mean(m14_1_df."Rho[1,2]")
	μ_a = mean(m14_1_df.a)
	μ_b = mean(m14_1_df.b)
	σ₁ = mean(m14_1_df."σ_cafe[1]")
	σ₂ = mean(m14_1_df."σ_cafe[2]")
	
	# draw ellipses
	dt = acos(ρ*σ₁*σ₂)
	chi = Chisq(2)
	
	for l ∈ (0.1, 0.3, 0.5, 0.8, 0.99)
	    scale = sqrt(quantile(chi, l))
	    xₜ(t) = scale*σ₁^2*cos(t + dt/2) + μ_a
	    yₜ(t) = scale*σ₂^2*cos(t - dt/2) + μ_b
	
	    plot!(xₜ, yₜ, 0, 2π, c=:black, alpha=0.3)
	end
	
	p_cafe
end

# ╔═╡ 5188efaa-8a5e-4fac-9e7f-bf048079ab9c
md" #### Code 14.16"

# ╔═╡ 23810d4f-1dd6-4c16-bf71-231b8b4b726d
let
	wait_morning_1 = a1
	wait_afternoon_1 = a1 .+ b1
	wait_morning_2 = a2
	wait_afternoon_2 = a2 .+ b2
	
	global p = scatter(wait_morning_1, wait_afternoon_1, xlab="morning wait", ylab="afternoon wait")
	scatter!(wait_morning_2, wait_afternoon_2, mc=:white)
	
	plot!(x -> x, s=:dash)
	
	for (x1,y1,x2,y2) ∈ zip(wait_morning_1, wait_afternoon_1, wait_morning_2, wait_afternoon_2)
	    plot!([x1,x2], [y1,y2], c=:black)
	end
	p
end

# ╔═╡ 91780708-0829-4beb-8144-7fc5b7e76b3f
md" #### Code 14.17"

# ╔═╡ cdeba577-a184-4bb5-97a9-6acf3ff516a3
let
	Random.seed!(1)
	
	# posterior mean
	ρ = mean(m14_1_df."Rho[1,2]")
	μ_a = mean(m14_1_df.a)
	μ_b = mean(m14_1_df.b)
	σ₁ = mean(m14_1_df."σ_cafe[1]")
	σ₂ = mean(m14_1_df."σ_cafe[2]")

	Σ = [[σ₁^2, σ₁*σ₂*ρ] [σ₁*σ₂*ρ, σ₂^2]]
	μ = [μ_a, μ_b]
	v = rand(MvNormal(μ, Σ), 10^4)
	v[2,:] += v[1,:]
	Σ₂ = cov(v')
	μ₂ = [μ_a, μ_a+μ_b]
	
	# draw ellipses
	dt = acos(Σ₂[1,2])
	chi = Chisq(2)
	
	for l ∈ (0.1, 0.3, 0.5, 0.8, 0.99)
	    scale = sqrt(quantile(chi, l))
	    xₜ(t) = scale*Σ₂[1,1]*cos(t + dt/2) + μ₂[1]
	    yₜ(t) = scale*Σ₂[2,2]*cos(t - dt/2) + μ₂[2]
	
	    plot!(xₜ, yₜ, 0, 2π, c=:black, alpha=0.3)
	end
	
	p
end

# ╔═╡ e049f64c-02e0-4087-b341-eaf437c5ba49
md" ## 14.2 Advanced varying slopes."

# ╔═╡ 687aec56-682c-4153-ac35-58926c03ea26
md" #### Code 14.18"

# ╔═╡ 0c918f86-e896-4b19-8cf8-a77203c3a61b
begin
	chimpanzees = CSV.read(sr_datadir("chimpanzees.csv"), DataFrame)
	chimpanzees.treatment = 1 .+ chimpanzees.prosoc_left .+ 2*chimpanzees.condition;
	chimpanzees.block_id = chimpanzees.block;
end;

# ╔═╡ fee0f70c-9aa5-4b60-a285-6923fedbd06c
@model function m14_2(L, tid, actor, block_id)
    tid_len = length(levels(tid))
    act_len = length(levels(actor))
    blk_len = length(levels(block_id))
    g ~ filldist(Normal(), tid_len)

    σ_actor ~ filldist(Exponential(), tid_len)
    ρ_actor ~ LKJ(tid_len, 2)
    Σ_actor = (σ_actor .* σ_actor') .* ρ_actor
    alpha ~ filldist(MvNormal(zeros(tid_len), Σ_actor), act_len)

    σ_block ~ filldist(Exponential(), tid_len)
    ρ_block ~ LKJ(tid_len, 2)
    Σ_block = (σ_block .* σ_block') .* ρ_block
    beta ~ filldist(MvNormal(zeros(tid_len), Σ_block), blk_len)
    
    for i ∈ eachindex(L)
        p = logistic(g[tid[i]] + alpha[tid[i], actor[i]] + beta[tid[i], block_id[i]])
        L[i] ~ Bernoulli(p)
    end
end

# ╔═╡ 57050894-2fcb-48e0-bd41-14502da1b5ae
begin
	Random.seed!(123)
	m14_2_ch = sample(m14_2(chimpanzees.pulled_left, chimpanzees.treatment, chimpanzees.actor, chimpanzees.block_id),
		HMC(0.01, 10), 1000);
	m14_2_df = DataFrame(m14_2_ch);
end

# ╔═╡ 2a1dd0cf-5a09-43bd-b7bf-5c93bea3bf3b
md" #### Code 14.19"

# ╔═╡ d15e5041-6021-483e-b059-7170077c2ea6
@model function m14_3(L, tid, actor, block_id)
    tid_len = length(levels(tid))
    act_len = length(levels(actor))
    blk_len = length(levels(block_id))
    g ~ filldist(Normal(), tid_len)

    σ_actor ~ filldist(Exponential(), tid_len)
    # LKJCholesky is not usable in Turing: https://github.com/TuringLang/Turing.jl/issues/1629
    ρ_actor ~ LKJ(tid_len, 2)    
    ρ_actor_L = cholesky(Symmetric(ρ_actor)).L
    z_actor ~ filldist(MvNormal(zeros(tid_len), 1), act_len)
    alpha = (σ_actor .* ρ_actor_L) * z_actor

    σ_block ~ filldist(Exponential(), tid_len)
    ρ_block ~ LKJ(tid_len, 2)
    ρ_block_L = cholesky(Symmetric(ρ_block)).L        
    z_block ~ filldist(MvNormal(zeros(tid_len), 1), blk_len)
    beta = (σ_block .* ρ_block_L) * z_block

    for i ∈ eachindex(L)
        p = logistic(g[tid[i]] + alpha[tid[i], actor[i]] + beta[tid[i], block_id[i]])
        L[i] ~ Bernoulli(p)
    end
end

# ╔═╡ 425eea72-221d-4990-90f2-a6654f4e98e1
md"

!!! note

Hm, this is less stable and slower than m14_2...
So if you know how to improve it - PR or open the issue in the repo"

# ╔═╡ ccc7d18c-08ed-41a7-a556-26030a19d6c5
begin
	Random.seed!(123)
	m14_3_ch = sample(m14_3(chimpanzees.pulled_left, chimpanzees.treatment, chimpanzees.actor, chimpanzees.block_id), 
    	HMC(0.01, 10), 1000)
	CHNS(m14_2_ch)
end

# ╔═╡ 6b73a8a5-c542-4a10-8cf1-0a7e40c3cc98
md" #### Code 14.20"

# ╔═╡ bc41a9ab-c960-4d17-9782-f01721337aea
let
	t = DataFrame(ess_rhat(m14_2_ch));
	ess_2 = t.ess[occursin.(r"σ|ρ", string.(t.parameters))]
	ess_2 = filter(v -> !isnan(v), ess_2);
	t = DataFrame(ess_rhat(m14_3_ch));
	ess_3 = t.ess[occursin.(r"σ|ρ", string.(t.parameters))]
	ess_3 = filter(v -> !isnan(v), ess_3);
	
	bounds = extrema([ess_2; ess_3]) .+ (-20, 20)
	scatter(ess_2, ess_3, xlim=bounds, ylim=bounds, xlab="centered (default)", ylab="non-centered (cholesky)")
	plot!(identity)
end

# ╔═╡ 40ec9daa-913e-4b2c-a0aa-3e70a2de7c71
md" #### Code 14.21"

# ╔═╡ b938e63c-29b1-4e49-90af-f9b2ab2e5bdd
begin
	m14_3_df = DataFrame(m14_3_ch)
	describe(m14_3_df[:, r"σ"])
end

# ╔═╡ e5428608-cfc2-4fb0-873d-9a4944fb27f6
md" #### Code 14.22"

# ╔═╡ 0d6baf1d-2e2a-4742-944a-072bbd24b628
md"

!!! note

The results for both models 2 and 3 are weird and mismatch with the book. So, something is wrong here.
Put below both link functions for experimentations.

Plot is from model 2, because model 3 is totally off."

# ╔═╡ ef730308-7bf8-489a-9ab5-a869efe6083a
let
	gd = groupby(chimpanzees, [:actor, :treatment])
	c = combine(gd, :pulled_left => mean => :val)
	global pl = unstack(c, :actor, :treatment, :val);
end

# ╔═╡ 1c4b4cb8-12c0-4fb3-8523-062c7028f746
l_fun = (r, (ai, ti)) -> begin
    bi = 5
    g = get(r, "g[$ti]", missing)
    
    σ_actor = get(r, "σ_actor[$ti]", missing)
    ρ_actor = reshape(collect(r[r"ρ_actor"]), (4, 4))
    ρ_actor_L = cholesky(Symmetric(ρ_actor)).L
    z_actor = reshape(collect(r[r"z_actor"]), (4, 7))
    alpha = (σ_actor .* ρ_actor_L) * z_actor
    a = alpha[ti, ai]
    
    σ_block = get(r, "σ_block[$ti]", missing)
    ρ_block = reshape(collect(r[r"ρ_block"]), (4, 4))
    ρ_block_L = cholesky(Symmetric(ρ_block)).L
    z_block = reshape(collect(r[r"z_block"]), (4, 6))
    beta = (σ_block .* ρ_block_L) * z_block
    b = beta[ti, bi]
    
    logistic(g + a + b)
end

# ╔═╡ 9610477f-8311-4cbf-8544-b9f4836926d6
#p_post = link(m14_3_df, l_fun, Iterators.product(1:7, 1:4))

l_fun2 = (r, (ai, ti)) -> begin
    bi = 5
    g = get(r, "g[$ti]", missing)
    a = get(r, "alpha[$ti,$ai]", missing)
    b = get(r, "beta[$ti,$bi]", missing)
    logistic(g + a + b)
end

# ╔═╡ e5216885-c156-4454-ac96-eed454435695
p_post = link(m14_2_df, l_fun2, Iterators.product(1:7, 1:4))

# ╔═╡ 6e47526a-310a-4a3e-aab5-0fe71c6c7814
let
	p_μ = map(mean, p_post)
	p_ci = map(PI, p_post);
	rel_ci = map(idx -> (p_μ[idx]-p_ci[idx][1], p_ci[idx][2]-p_μ[idx]), CartesianIndices(p_ci));
	
	n_names = ["R/N", "L/N", "R/P", "L/P"]
	p = plot(ylims=(0, 1.1), ylab="proportion left lever", showaxis=:y, xticks=false)
	hline!([0.5], c=:gray, s=:dash)
	
	# raw data
	for actor in 1:7
	    ofs = (actor-1)*4
	    actor > 1 && vline!([ofs+0.5], c=:gray)
	    plot!([ofs+1,ofs+3], collect(pl[actor,["1","3"]]), lw=2, m=:o, c=:black)
	    plot!([ofs+2,ofs+4], collect(pl[actor,["2","4"]]), lw=2, m=:o, c=:black)
	    anns = [
	        (ofs+idx, pl[actor,string(idx)]+.04, (name, 8))
	        for (idx,name) ∈ enumerate(n_names)
	    ]
	    actor != 2 && annotate!(anns)
	end
	
	annotate!([
	    (2.5 + (idx-1)*4, 1.1, ("actor $idx", 8))
	    for idx ∈ 1:7
	])
	
	# posterior predictions
	for actor in 1:7
	    ofs = (actor-1)*4
	    actor > 1 && vline!([ofs+0.5], c=:gray)
	    err = [rel_ci[actor,1], rel_ci[actor,3]]
	    plot!([ofs+1,ofs+3], collect(p_μ[actor,[1,3]]), err=err, lw=2, m=:o, c=:blue)
	    err = [rel_ci[actor,2], rel_ci[actor,4]]    
	    plot!([ofs+2,ofs+4], collect(p_μ[actor,[2,4]]), err=err, lw=2, m=:o, c=:blue)
	end
	
	p
end

# ╔═╡ 944c59b8-ecf6-4ec3-9e81-7ffb8692db36
md" ## 14.3 Instruments and causal designs."

# ╔═╡ 233cd90e-463d-49fb-9a54-34e2750c2a8c
md" #### Code 14.23"

# ╔═╡ b28587b9-5e60-42b0-9874-d47b3e03b97c
let
	Random.seed!(73)
	N = 500
	U_sim = rand(Normal(), N)
	Q_sim = rand(1:4, N)
	E_sim = [rand(Normal(μ)) for μ ∈ U_sim .+ Q_sim]
	W_sim = [rand(Normal(μ)) for μ ∈ U_sim .+ 0*Q_sim]
	
	global dat_sim1 = DataFrame(
	    W=standardize(ZScoreTransform, W_sim),
	    E=standardize(ZScoreTransform, E_sim),
	    Q=standardize(ZScoreTransform, float.(Q_sim)),
	)
end

# ╔═╡ 4fb1f19d-2802-4a0f-bb35-1e21c0be68f3
md" #### Code 14.24"

# ╔═╡ aad09fb4-1fee-454f-bdad-65b2a93255ca
@model function m14_4(W, E)
    σ ~ Exponential()
    aW ~ Normal(0, 0.2)
    bEW ~ Normal(0, 0.5)
    μ = @. aW + bEW * E
    W ~ MvNormal(μ, σ)
end

# ╔═╡ 0d451930-833e-4050-b2f9-d5f997fb0098
begin
	m14_4_ch = sample(m14_4(dat_sim1.W, dat_sim1.E),
	    NUTS(), 1000)
	m14_4_df = DataFrame(m14_4_ch)
	describe(m14_4_df)
end

# ╔═╡ 3878adef-0d96-4e73-ab37-4e761e01e0ac
md" #### Code 14.25"

# ╔═╡ 500fb09f-f952-46a0-8d63-98279e630882
@model function m14_5(W, E, Q)
    σ ~ Exponential()
    aW ~ Normal(0, 0.2)
    bEW ~ Normal(0, 0.5)
    bQW ~ Normal(0, 0.5)
    μ = @. aW + bEW * E + bQW * Q
    W ~ MvNormal(μ, σ)
end

# ╔═╡ 95d072b9-5940-49ce-97a8-8dde13c25939
begin
	m14_5_ch = sample(m14_5(dat_sim1.W, dat_sim1.E, dat_sim1.Q),
	    NUTS(), 1000)
	m14_5_df = DataFrame(m14_5_ch)
	describe(m14_5_df)
end

# ╔═╡ 1f3d5631-8261-4c3b-8e34-761d23eb5e2e
md" #### Code 14.26"

# ╔═╡ a856d989-5fb6-4148-ac7d-2ce52952c1a3
@model function m14_6(W, E, Q, WE)
    σ ~ filldist(Exponential(), 2)
    ρ ~ LKJ(2, 2)
    aW ~ Normal(0, 0.2)
    aE ~ Normal(0, 0.2)
    bEW ~ Normal(0, 0.5)
    bQE ~ Normal(0, 0.5)
    μW = @. aW + bEW*E
    μE = @. aW + bQE*Q
    Σ = (σ .* σ') .* ρ
    for i ∈ eachindex(WE)
        WE[i] ~ MvNormal([μW[i], μE[i]], Σ)
    end
end

# ╔═╡ 9a7ebb96-5ce1-4093-8c85-9522edebab0e
begin
	Random.seed!(1)
	# need to combine W and E here (Turing vars limitation)
	WE = [[w,e] for (w,e) ∈ zip(dat_sim1.W, dat_sim1.E)]
	m14_6_ch = sample(m14_6(dat_sim1.W, dat_sim1.E, dat_sim1.Q, WE), NUTS(200, 0.65, init_ϵ=0.003), 1000)
	m14_6_df = DataFrame(m14_6_ch);
end

# ╔═╡ 48d158ed-a25f-4ebe-8863-6c92071c251c
let
	# Drop cols with zero variance
	df = m14_6_df[:, Not("ρ[1,1]")][:, Not("ρ[2,2]")]
	describe(df)
end

# ╔═╡ 199ae6d5-7306-4802-bd18-4232fdb98f5f
md" #### Code 14.28"

# ╔═╡ aee19854-5c3d-4a9e-8639-91843800d1b6
let
	Random.seed!(73)
	
	N = 500
	U_sim = rand(Normal(), N)
	Q_sim = rand(1:4, N)
	E_sim = [rand(Normal(μ)) for μ ∈ U_sim .+ Q_sim]
	W_sim = [rand(Normal(μ)) for μ ∈ -U_sim .+ 0.2*Q_sim]
	
	global dat_sim2 = DataFrame(
	    W=standardize(ZScoreTransform, W_sim),
	    E=standardize(ZScoreTransform, E_sim),
	    Q=standardize(ZScoreTransform, float.(Q_sim)),
	);
end

# ╔═╡ e01fa7a9-16ac-498e-97d3-8c9fe1255f43
md" #### Code 14.29"

# ╔═╡ 48fe8dff-dc91-4996-bfe7-911d2a09b4c8
md"

!!! note

Not implemented in dagitty.jl yet."

# ╔═╡ aef2347f-dea0-4bbb-864c-c7d0d09814c6
let
	g = DAG(:Q => :E, :U => :E, :E => :W, :E => :W)
end

# ╔═╡ 011c2852-91b3-4745-b453-e0ffcdba7f23
md" ## 14.4 Social relations as correlated varying effects."

# ╔═╡ 6916ba32-8ebc-4aa7-adf7-f6717fde60e2
md" #### Code 14.30"

# ╔═╡ afaf9be1-a557-4b73-9279-f9cf12525958
begin
	kl_dyads = CSV.read(sr_datadir("KosterLeckie.csv"), DataFrame)
	describe(kl_dyads)
end

# ╔═╡ 1c03fd5e-0c6a-4118-b11a-0f4d2392f923
md" #### Code 14.31"

# ╔═╡ 1cd43a15-d361-4e4a-84e2-f1464659a654
# +
kl_data = (
    N = nrow(kl_dyads), 
    N_households = maximum(kl_dyads.hidB),
    did = kl_dyads.did,
    hidA = kl_dyads.hidA,
    hidB = kl_dyads.hidB,
    giftsAB = kl_dyads.giftsAB,
    giftsBA = kl_dyads.giftsBA,
)

# ╔═╡ 2be10e8e-9a44-4b0b-bd81-7e4cdc999afc
@model function m14_7(N, N_households, hidA, hidB, did, giftsAB, giftsBA)
    a ~ Normal()
    ρ_gr ~ LKJ(2, 4)
    σ_gr ~ filldist(Exponential(), 2)
    Σ = (σ_gr .* σ_gr') .* ρ_gr
    gr ~ filldist(MvNormal(Σ), N_households)
    
    # dyad effects (use 2 z values)
    z₁ ~ filldist(Normal(), N)
    z₂ ~ filldist(Normal(), N)
    z = [z₁ z₂]'
    σ_d ~ Exponential()
    ρ_d ~ LKJ(2, 8)
    L_ρ_d = cholesky(Symmetric(ρ_d)).L
    d = (σ_d .* L_ρ_d) * z

    λ_AB = exp.(a .+ gr[1, hidA] .+ gr[2, hidB] .+ d[1, did])
    λ_BA = exp.(a .+ gr[1, hidB] .+ gr[2, hidA] .+ d[2, did])
    for i ∈ eachindex(giftsAB)
        giftsAB[i] ~ Poisson(λ_AB[i])
        giftsBA[i] ~ Poisson(λ_BA[i])
    end
    return d
end

# ╔═╡ ee84ac32-58f1-49bc-b705-fbc709e0f11a
begin
	model = m14_7(
	    kl_data.N, kl_data.N_households, kl_data.hidA, kl_data.hidB, 
	    kl_data.did, kl_data.giftsAB, kl_data.giftsBA
	)
	m14_7_ch = sample(model, 
	    NUTS(1000, 0.65, init_ϵ=0.025), 
	    1000)
	m14_7_df = DataFrame(m14_7_ch);
end

# ╔═╡ 8448e4cd-0bce-4426-85db-72a40a279b1d
md" #### Code 14.32"

# ╔═╡ 640ac000-dd55-4dda-846f-303404360ab5
describe(m14_7_df[:, r"_gr\[(1,2|2,1|1|2)\]"])

# ╔═╡ 038bc011-279e-492a-bbfc-074630cdff3e
md" #### Code 14.33"

# ╔═╡ c8d7399e-ad0d-4835-8fa0-1c8ad66e7d34
let
	g = [
	    m14_7_df.a .+ m14_7_df[!,"gr[1,$i]"]
	    for i ∈ 1:25
	]
	r = [
	    m14_7_df.a .+ m14_7_df[!,"gr[2,$i]"]
	    for i ∈ 1:25
	]
	g = hcat(g...)'
	r = hcat(r...)';
	Eg_μ = mean(eachcol(exp.(g)))
	Er_μ = mean(eachcol(exp.(r)));
	
	# Code 14.34
	
	# +
	plot(xlim=(0, 8.6), ylim=(0,8.6), xlab="generalized giving", ylab="generalized receiving")
	plot!(x -> x, c=:black, s=:dash)
	
	for i ∈ 1:25
	    gi = exp.(g[i,:])
	    ri = exp.(r[i,:])
	    Σ = cov([gi ri])
	    μ = [mean(gi), mean(ri)]
	
	    dt = acos(Σ[1,2])
	    xₜ(t) = Σ[1,1]*cos(t + dt/2) + μ[1]
	    yₜ(t) = Σ[2,2]*cos(t - dt/2) + μ[2]
	
	    plot!(xₜ, yₜ, 0, 2π, c=:black, lw=1)
	end
	
	scatter!(Eg_μ, Er_μ, c=:white, msw=1.5)
end

# ╔═╡ 1512c0e0-bee0-4fce-9062-39c60e43f3f8
md" #### Code 14.35"

# ╔═╡ 93a5e590-0ab0-4a13-8944-226258424f7d
describe(m14_7_df[:, r"_d"])

# ╔═╡ aba1a47c-534b-4ce0-b41f-d68c3331207d
md" #### Code 14.36"

# ╔═╡ 2ff1b0f7-4a07-467f-bf0d-069c5b766ed1
#
# Illustrates `generated_quantities` trick to extract values returned from the model

let
	ch = Turing.MCMCChains.get_sections(m14_7_ch, :parameters)
	d_vals = generated_quantities(model, ch)
	
	d_y1 = [r[1,:] for r in d_vals]
	d_y1 = hcat(d_y1...)
	d_y1 = mean.(eachrow(d_y1))
	
	d_y2 = [r[2,:] for r in d_vals]
	d_y2 = hcat(d_y2...)
	d_y2 = mean.(eachrow(d_y2))
	
	scatter(d_y1, d_y2)
end

# ╔═╡ 0f77c94d-011e-4a55-99a1-0ee67acf63fb
md" ## 14.5 Continuous categories and the Gaussian process."

# ╔═╡ Cell order:
# ╠═95965ab2-f250-479f-9d5d-e478dba1fe35
# ╠═4458903f-4717-4603-a707-4d2fd587470b
# ╠═d577b8c8-3670-497a-a6fe-41208a8e8501
# ╠═90028541-7308-476f-b51b-850cd6ad39c4
# ╟─1217e4ab-8e97-4b73-9493-a58c7d8165fe
# ╟─100c3470-b7e4-4890-9347-851b78f5202f
# ╟─464ddfd5-d06b-4615-bd23-e63e03a29f1e
# ╠═fe81bd6a-c7bf-437d-8ca6-f21f9b45b41f
# ╟─b6c8f7af-8e39-4792-a62c-88a41876345d
# ╠═8ee2eec3-3214-48d1-8604-a54ba1a67baf
# ╟─573253a5-10d4-458e-9451-89c763658810
# ╠═fa9dff6f-3182-4093-bbfd-8deea404eab0
# ╟─2367f3cf-c568-4f6c-82c6-b7ae42408404
# ╠═8185285a-1640-45e9-8c80-b17b6c960804
# ╟─40634f53-8104-4791-877c-2d498cab9570
# ╠═885a0261-3541-4389-beb7-231c77e2dfcc
# ╟─6837fff8-0c29-4ba1-9b3b-b21c23d861eb
# ╠═1f30d601-4e00-44dd-a6ab-363c431ca875
# ╠═de06a783-c2ff-4f45-8a1e-02977ad31e26
# ╟─30be727e-86e6-41cf-be58-e0b1ffb372c2
# ╠═5ca07e4a-bf45-4e3c-9a5b-3b9b7af38020
# ╟─f95ec79a-0263-46e4-b82e-a274cb7ee7fb
# ╟─81f57cdb-36f8-43bb-86a5-19bac5cd88b1
# ╠═78fe51f6-beea-4155-905d-2d5e6e7d66d7
# ╟─4ea3a429-e457-48a1-a66a-ea66ff21d22c
# ╠═8bb5afba-3f59-4776-9a7b-8df9fed58a90
# ╟─5188efaa-8a5e-4fac-9e7f-bf048079ab9c
# ╠═23810d4f-1dd6-4c16-bf71-231b8b4b726d
# ╟─91780708-0829-4beb-8144-7fc5b7e76b3f
# ╠═cdeba577-a184-4bb5-97a9-6acf3ff516a3
# ╟─e049f64c-02e0-4087-b341-eaf437c5ba49
# ╟─687aec56-682c-4153-ac35-58926c03ea26
# ╠═0c918f86-e896-4b19-8cf8-a77203c3a61b
# ╠═fee0f70c-9aa5-4b60-a285-6923fedbd06c
# ╠═57050894-2fcb-48e0-bd41-14502da1b5ae
# ╟─2a1dd0cf-5a09-43bd-b7bf-5c93bea3bf3b
# ╠═d15e5041-6021-483e-b059-7170077c2ea6
# ╟─425eea72-221d-4990-90f2-a6654f4e98e1
# ╠═ccc7d18c-08ed-41a7-a556-26030a19d6c5
# ╟─6b73a8a5-c542-4a10-8cf1-0a7e40c3cc98
# ╠═bc41a9ab-c960-4d17-9782-f01721337aea
# ╟─40ec9daa-913e-4b2c-a0aa-3e70a2de7c71
# ╠═b938e63c-29b1-4e49-90af-f9b2ab2e5bdd
# ╟─e5428608-cfc2-4fb0-873d-9a4944fb27f6
# ╟─0d6baf1d-2e2a-4742-944a-072bbd24b628
# ╠═ef730308-7bf8-489a-9ab5-a869efe6083a
# ╠═1c4b4cb8-12c0-4fb3-8523-062c7028f746
# ╠═9610477f-8311-4cbf-8544-b9f4836926d6
# ╠═e5216885-c156-4454-ac96-eed454435695
# ╠═6e47526a-310a-4a3e-aab5-0fe71c6c7814
# ╟─944c59b8-ecf6-4ec3-9e81-7ffb8692db36
# ╟─233cd90e-463d-49fb-9a54-34e2750c2a8c
# ╠═b28587b9-5e60-42b0-9874-d47b3e03b97c
# ╟─4fb1f19d-2802-4a0f-bb35-1e21c0be68f3
# ╠═aad09fb4-1fee-454f-bdad-65b2a93255ca
# ╠═0d451930-833e-4050-b2f9-d5f997fb0098
# ╟─3878adef-0d96-4e73-ab37-4e761e01e0ac
# ╠═500fb09f-f952-46a0-8d63-98279e630882
# ╠═95d072b9-5940-49ce-97a8-8dde13c25939
# ╟─1f3d5631-8261-4c3b-8e34-761d23eb5e2e
# ╠═a856d989-5fb6-4148-ac7d-2ce52952c1a3
# ╠═9a7ebb96-5ce1-4093-8c85-9522edebab0e
# ╠═48d158ed-a25f-4ebe-8863-6c92071c251c
# ╟─199ae6d5-7306-4802-bd18-4232fdb98f5f
# ╠═aee19854-5c3d-4a9e-8639-91843800d1b6
# ╟─e01fa7a9-16ac-498e-97d3-8c9fe1255f43
# ╟─48fe8dff-dc91-4996-bfe7-911d2a09b4c8
# ╠═aef2347f-dea0-4bbb-864c-c7d0d09814c6
# ╟─011c2852-91b3-4745-b453-e0ffcdba7f23
# ╟─6916ba32-8ebc-4aa7-adf7-f6717fde60e2
# ╠═afaf9be1-a557-4b73-9279-f9cf12525958
# ╟─1c03fd5e-0c6a-4118-b11a-0f4d2392f923
# ╠═1cd43a15-d361-4e4a-84e2-f1464659a654
# ╠═2be10e8e-9a44-4b0b-bd81-7e4cdc999afc
# ╠═ee84ac32-58f1-49bc-b705-fbc709e0f11a
# ╠═8448e4cd-0bce-4426-85db-72a40a279b1d
# ╠═640ac000-dd55-4dda-846f-303404360ab5
# ╠═038bc011-279e-492a-bbfc-074630cdff3e
# ╠═c8d7399e-ad0d-4835-8fa0-1c8ad66e7d34
# ╠═1512c0e0-bee0-4fce-9062-39c60e43f3f8
# ╠═93a5e590-0ab0-4a13-8944-226258424f7d
# ╟─aba1a47c-534b-4ce0-b41f-d68c3331207d
# ╠═2ff1b0f7-4a07-467f-bf0d-069c5b766ed1
# ╟─0f77c94d-011e-4a55-99a1-0ee67acf63fb
