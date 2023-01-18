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

# ╔═╡ 3ae8474f-90d0-4653-a336-0803f7705462
begin
	using GLM
	using CSV
	using Random
	using StatsBase
	using DataFrames
	using Dagitty
	using Turing
	using StatsPlots
	using StatisticalRethinking
	using StatisticalRethinkingPlots
	using Logging
end

# ╔═╡ 579065ea-bc4c-49a3-a3d9-673e49071dfe
md" ## The haunted DAG and the causal terror."

# ╔═╡ 648d7e4e-5575-4867-814f-2972fe6c6729
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

# ╔═╡ 8da6e45a-9d78-479a-a8bb-57641e9a7375
let
	default(labels=false)
	Logging.disable_logging(Logging.Warn)
end;

# ╔═╡ 299eb270-63f8-496f-a235-400f3859babd
md" #### Code 6.1"

# ╔═╡ 40d57102-ed75-4da4-bdb1-cc18243ea570
begin
	Random.seed!(1917)
	N = 200   # grant proposals
	p = 0.1   # proportion to select
	
	# uncorrelated newsworthiness and trustworthiness
	nw = rand(Normal(), N)
	tw = rand(Normal(), N)
	
	# select top 10% of combined score
	s = nw .+ tw
	q = quantile(s, 1-p)
	selected = s .>= q
	cor(tw[selected], nw[selected])
end

# ╔═╡ 399bb4a5-23f2-4c9e-bae0-61ae8c668bed
let
	scatter(nw[.!selected], tw[.!selected]; xlab="newsworthiness", ylab="trustworthiness", label="not selected")
	scatter!(nw[selected], tw[selected]; label="selected")
end

# ╔═╡ 008c3130-6d49-45dd-a1da-c125c624b23b
md" ## 6.1 Multicollinearity."

# ╔═╡ 54c70b49-beee-4ac3-92f0-d7cd2816a6be
md" #### Code 6.2"

# ╔═╡ 0ef7a233-2a4e-4bc1-bd46-7387165b5be7
let
	Random.seed!(100)
	N = 100
	height = rand(Normal(10, 2), N)
	leg_prop = rand(Uniform(0.4, 0.5), N)
	leg_left = leg_prop .* height .+ rand(Normal(0, 0.02), N)
	leg_right = leg_prop .* height .+ rand(Normal(0, 0.02), N)
	global height_df = DataFrame(:height => height, :leg_left => leg_left, :leg_right => 		leg_right);
end

# ╔═╡ 63369b75-a94c-4365-8171-79c72051eaa1
md" #### Code 6.3"

# ╔═╡ 318c8900-a52d-4338-98d0-9dc2f1073cea
@model function model_m6_1(leg_left, leg_right, height)
    a ~ Normal(10, 100)
    bl ~ Normal(2, 10)
    br ~ Normal(2, 10)
    μ = @. a + bl * leg_left + br * leg_right
    σ ~ Exponential(1)
    height ~ MvNormal(μ, σ)
end

# ╔═╡ 40e28711-0a2e-40d2-a584-ae22f4de621a
begin
	m6_1 = sample(model_m6_1(height_df.leg_left, height_df.leg_right, height_df.height), NUTS(), 1000)
	m6_1_df = DataFrame(m6_1)
	PRECIS(m6_1_df)
end

# ╔═╡ 999c880d-596e-4ad2-ab4e-f2cc364ac99c
md" #### Code 6.4"

# ╔═╡ d9bd103f-2cc5-48cc-89e0-99c08d2897bf
coeftab_plot(m6_1_df)

# ╔═╡ aec96c68-fa64-4548-b68f-b6b0f2406182
md" #### Code 6.5"

# ╔═╡ f4eecae9-0374-4845-b673-c11bdedc2819
scatter(m6_1_df.br, m6_1_df.bl; alpha=0.1)

# ╔═╡ c7988aad-113e-4df4-95ed-90acfbaf60a2
md" #### Code 6.6"

# ╔═╡ bda99113-bd50-499f-aa7c-e6eea69d48da
@df m6_1_df density(:br + :bl; lw=2, xlab="sum of bl and br")

# ╔═╡ ea3ecd44-ccfd-41d8-91b9-10301a49523f
md" #### Code 6.7"

# ╔═╡ b02ec264-69ad-4d2c-a847-a9db6396631c
@model function model_m6_2(leg_left, height)
    a ~ Normal(10, 100)
    bl ~ Normal(2, 10)
    μ = @. a + bl * leg_left
    σ ~ Exponential(1)
    height ~ MvNormal(μ, σ)
end

# ╔═╡ 04be25da-56cd-49c9-ae84-466f6fe672d8
begin
	m6_2 = sample(model_m6_2(height_df.leg_left, height_df.height), NUTS(), 1000)
	m6_2_df = DataFrame(m6_2)
	PRECIS(m6_2_df)
end

# ╔═╡ 5a81fd3c-fc60-47a5-9da9-2e40eeb30251
std(m6_1_df.bl), std(m6_1_df.br), std(m6_1_df.bl + m6_1_df.br)

# ╔═╡ 0e9cd9a9-89b3-465a-8936-f10ef02eb4d0
md" #### Code 6.8"

# ╔═╡ 4e8ea51e-dedb-4f21-8661-4220b6a6227b
begin
	milk = DataFrame(CSV.File(sr_datadir("milk.csv"),  missingstring="NA"))
	
	# get rid of dots in column names
	rename!(n -> replace(n, "." => "_"), milk)
	
	milk[!,:K] = standardize(ZScoreTransform, milk.kcal_per_g)
	milk[!,:F] = standardize(ZScoreTransform, milk.perc_fat)
	milk[!,:L] = standardize(ZScoreTransform, milk.perc_lactose);
end;

# ╔═╡ abfee492-d400-4543-a7a4-a3382bff4087
md" #### Code 6.9"

# ╔═╡ a9b652e5-3376-43c0-b962-d905810c82e8
@model function model_m6_3(F, K)
    a ~ Normal(0, 0.2)
    bF ~ Normal(0, 0.5)
    μ = @. a + F * bF
    σ ~ Exponential(1)
    K ~ MvNormal(μ, σ)
end

# ╔═╡ 020f1bbd-3f89-4256-90c0-b595712128bc
begin
	m6_3 = sample(model_m6_3(milk.F, milk.K), NUTS(), 1000)
	m6_3_df = DataFrame(m6_3)
end

# ╔═╡ d17704eb-1212-42a3-bcc4-f59bf3231333
@model function model_m6_4(L, K)
    a ~ Normal(0, 0.2)
    bL ~ Normal(0, 0.5)
    μ = @. a + L * bL
    σ ~ Exponential(1)
    K ~ MvNormal(μ, σ)
end

# ╔═╡ badb90cf-aaa4-4e69-92da-54ddf020b473
begin
	m6_4 = sample(model_m6_4(milk.L, milk.K), NUTS(), 1000)
	m6_4_df = DataFrame(m6_4)
end

# ╔═╡ 1755c963-e12b-4f5a-b670-51e7cf3257f9
PRECIS(m6_3_df)

# ╔═╡ fcc55095-6ae2-40eb-81bd-639a41550f0a
PRECIS(m6_4_df)

# ╔═╡ 61e1e17d-1347-44d6-95b9-02ecd9f58450
md" #### Code 6.10"

# ╔═╡ b6208f9a-b80a-4b0a-8833-7a3154b0b96c
@model function model_m6_5(F, L, K)
    a ~ Normal(0, 0.2)
    bF ~ Normal(0, 0.5)
    bL ~ Normal(0, 0.5)
    μ = @. a + F * bF + L * bL
    σ ~ Exponential(1)
    K ~ MvNormal(μ, σ)
end

# ╔═╡ cf869b8f-2c35-4ba3-9546-2d15cbdaf6e0
begin
	m6_5 = sample(model_m6_5(milk.F, milk.L, milk.K), NUTS(), 1000)
	m6_5_df = DataFrame(m6_5)
	PRECIS(m6_5_df)
end

# ╔═╡ 8e1f339c-cf39-4b73-b0aa-fbd37c73e81d
md" #### Code 6.11"

# ╔═╡ 932dc19e-080f-451e-88a7-653f2ba82b1d
@df milk corrplot([:kcal_per_g :perc_fat :perc_lactose]; seriestype=:scatter, bins=10, grid=false)

# ╔═╡ a2b8eeb5-67df-4023-8631-6b001feb20c7
md" #### Code 6.12"

# ╔═╡ 7d79c5ed-e6b6-4635-9cfa-8847bd2da8ae
# get mean stderr for linear model's scale
function stderr_for_r(r)
    σ = sqrt(1-r^2)*var(milk.perc_fat)
    fat_scaled = r .* milk.perc_fat
    stderr_x = [
        begin
            x = milk.perc_fat .+ rand(MvNormal(fat_scaled, σ))
            # add the intercept to the model
            X = hcat(ones(length(x)), x)
            m = lm(X, milk.kcal_per_g)
            stderror(m)[2]
        end
        for _ in 1:100
    ]
    s = mean(stderr_x)
end

# ╔═╡ fc0f6c2a-e625-4fa6-bcf7-1bc64f617c27
let
	r_seq = range(0, 0.99; step=0.01)
	s = stderr_for_r.(r_seq)
	plot(r_seq, s; lw=2, xlab="correlation")
end

# ╔═╡ 2dc14c15-86b2-4b78-965a-e63fac38db6f
md" ## 6.2 Post-treatment bias"

# ╔═╡ 03e1487b-321f-41a8-b2e8-a27d53b29fa6
md" #### Code 6.13"

# ╔═╡ b89b0a75-3dc3-443e-af41-9b17df5dd649
let
	Random.seed!(70)
	# number of plants
	N = 100
	h0 = rand(Normal(10, 2), N)
	treatment = repeat(0:1, inner=div(N, 2))
	fungus = [rand(Binomial(1, 0.5 - treat*0.4)) for treat in treatment]
	h1 = h0 .+ rand(MvNormal(5 .- 3 .* fungus, 1))
	global fungus_df = DataFrame(:h0 => h0, :h1 => h1, :treatment => treatment, :fungus => fungus)
	PRECIS(fungus_df)
end

# ╔═╡ 83264c24-cdbf-44a7-b194-6f4380ce3a12
md" #### Code 6.14"

# ╔═╡ 2e9c86d6-33af-4fed-ba76-304b5e8ba6fe
let
	sim_p = rand(LogNormal(0, 0.25), 10_000)
	PRECIS(DataFrame(:sim_p => sim_p))
end

# ╔═╡ ab057169-5969-4c80-8028-61dabf9b1f1a
md" #### Code 6.15"

# ╔═╡ 4005549c-b251-4485-9107-c5bfee420bdf
@model function model_m6_6(h0, h1)
    p ~ LogNormal(0, 0.25)
    σ ~ Exponential(1)
    μ = h0 .* p
    h1 ~ MvNormal(μ, σ)
end

# ╔═╡ ef69951e-cfdc-487b-94d6-4b21217b4a8c
begin
	m6_6 = sample(model_m6_6(fungus_df.h0, fungus_df.h1), NUTS(), 1000)
	m6_6_df = DataFrame(m6_6)
	PRECIS(m6_6_df)
end

# ╔═╡ 52f0e6aa-f55e-493c-b8a1-f681170b5ef8
md" #### Code 6.16"

# ╔═╡ f8b2e102-df8a-4198-a242-58ab278d6ccc
@model function model_m6_7(h0, treatment, fungus, h1)
    a ~ LogNormal(0, 0.2)
    bt ~ Normal(0, 0.5)
    bf ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    p = @. a + bt*treatment + bf*fungus
    μ = h0 .* p
    h1 ~ MvNormal(μ, σ)
end

# ╔═╡ 43526e84-61d1-40ae-8146-416b9406e7e4
begin
	m6_7 = sample(model_m6_7(fungus_df.h0, fungus_df.treatment, fungus_df.fungus, fungus_df.h1), NUTS(), 1000)
	m6_7_df = DataFrame(m6_7)
	PRECIS(m6_7_df)
end

# ╔═╡ fc43c304-9e82-4e4f-b56a-44cd7f664520
md" #### Code 6.17"

# ╔═╡ b507451b-b0e9-4786-80ab-fdae9c86f039
@model function model_m6_8(h0, treatment, h1)
    a ~ LogNormal(0, 0.2)
    bt ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    p = @. a + bt*treatment
    μ = h0 .* p
    h1 ~ MvNormal(μ, σ)
end

# ╔═╡ 7faa03d8-23b4-4792-9b08-080f357be27a
begin
	m6_8_2 = sample(model_m6_8(fungus_df.h0, fungus_df.treatment, fungus_df.h1), NUTS(), 1000)
	m6_8_df = DataFrame(m6_8_2)
	PRECIS(m6_8_df)
end

# ╔═╡ 08f49590-0dd0-40d4-ae07-74b90d923811
md" #### Code 6.18"

# ╔═╡ 052f4881-2db3-49d9-817f-c09c7abd7b80
let
	plant_dag = Dagitty.DAG(:H₀ => :H₁, :F => :H₁, :T => :F)
	drawdag(plant_dag, [2, 0, 1, 3], [0, 0, 0, 0])
end

# ╔═╡ 13242e24-b521-4181-9d57-d9dabe94c322
md" #### Code 6.19"

# ╔═╡ b79ebc48-f738-491e-8fe2-822c1828ad90
implied_conditional_independencies_min(plant_dag)

# ╔═╡ 8902c6c3-f99a-48b5-bbd9-5ce65dda504f
md" #### Code 6.20"

# ╔═╡ 695028d4-42a9-469b-94ae-0456268e4c62
let
	Random.seed!(70)
	# number of plants
	N = 1000
	h0 = rand(Normal(10, 2), N)
	treatment = repeat(0:1, inner=div(N, 2))
	M = rand(Bernoulli(), N)
	fungus = [
	    rand(Binomial(1, 0.5 - treat*0.4 + 0.4 * m)) 
	    for (treat, m) ∈ zip(treatment, M)
	]
	h1 = h0 .+ rand(MvNormal(5 .+ 3 .* M, 1))
	
	global fungus_df2 = DataFrame(:h0 => h0, :h1 => h1, :treatment => treatment, :fungus => fungus)
	PRECIS(fungus_df2)
end

# ╔═╡ 0ce498fd-cec5-4028-872c-88f8941a0e87
begin
	m6_7_2 = sample(model_m6_7(fungus_df2.h0, fungus_df2.treatment, fungus_df2.fungus, fungus_df2.h1), NUTS(), 1000)
	PRECIS(DataFrame(m6_7_2))
end

# ╔═╡ fe3d9e94-1ac1-48cf-9a58-ec409324f7f4
begin
	m6_8 = sample(model_m6_8(fungus_df2.h0, fungus_df2.treatment, fungus_df2.h1), NUTS(), 1000)
	PRECIS(DataFrame(m6_8))
end

# ╔═╡ 3659ca0e-5ce4-4227-a6e1-b4f250e6ba30
md" ## 6.3 Collider bias"

# ╔═╡ be9f02e7-3704-4764-bc1f-0adef06f2836
md" #### Code 6.21"

# ╔═╡ d0c867a1-e69e-4034-9dbb-80f87dc33a2a
begin
	happiness = sim_happiness(seed=1977, n_years=1000)
	precis(happiness)
end

# ╔═╡ a0194048-1517-4dcf-84af-dafbc0a0222c
let
	d_m = happiness[happiness.married .== 1,[:age,:happiness]]
	d_u = happiness[happiness.married .== 0,[:age,:happiness]]
	
	scatter(d_m.age, d_m.happiness; label="married", xlab="age", ylab="happiness")
	scatter!(d_u.age, d_u.happiness; c=:white)
end

# ╔═╡ e1269e2b-0f88-4137-b182-2f0c8f79f6b3
md" #### Code 6.22"

# ╔═╡ cbef9b36-bf1d-40c9-83d0-72879c683010
begin
	happiness2 = happiness[happiness.age .> 17,:]
	happiness2[!,:A] = @. (happiness2.age - 18) / (65-18);
end;

# ╔═╡ 0ad131a0-212a-4169-9a69-033d02cba337
md" #### Code 6.23"

# ╔═╡ 576750fa-4ff7-477b-8a17-2707c0b9b4c1
happiness2[!,:mid] = happiness2.married .+ 1;

# ╔═╡ a1832643-71b9-4728-ad53-cca2fae2da07
@model function model_m6_9(mid, A, happiness)
    a ~ MvNormal([0, 0], 1)
    bA ~ Normal(0, 2)
    μ = a[mid] .+ bA .* A
    σ ~ Exponential(1)
    happiness ~ MvNormal(μ, σ)
end

# ╔═╡ b8dadb5d-bc7a-4f35-8f7e-8eb803cc6c46
begin
	m6_9 = sample(model_m6_9(happiness2.mid, happiness2.A, happiness2.happiness), NUTS(), 1000)
	m6_9_df = DataFrame(m6_9)
	PRECIS(m6_9_df)
end

# ╔═╡ f676cc8f-8ec3-4e30-9a6d-493b65b894be
md" #### Code 6.24"

# ╔═╡ 791d6baf-d0bb-4c8f-997c-60e9d5e73e51
@model function model_m6_10(A, happiness)
    a ~ Normal()
    bA ~ Normal(0, 2)
    μ = a .+ bA .* A
    σ ~ Exponential(1)
    happiness ~ MvNormal(μ, σ)
end

# ╔═╡ a33bd041-06d5-4cd4-a88c-38bc8ff4c0cc
begin
	m6_10 = sample(model_m6_10(happiness2.A, happiness2.happiness), NUTS(), 1000)
	m6_10_df = DataFrame(m6_10)
	PRECIS(m6_10_df)
end

# ╔═╡ 8836c691-492a-4f76-9193-05389b69ba34
md" #### Code 6.25"

# ╔═╡ c81ca3a3-8ea5-4ccb-a420-75bb0970d5b3
begin
	b_GP = 1
	b_GC = 0
	b_PC = 1
	b_U = 2;
end;

# ╔═╡ b43cbceb-79be-4f66-9a4b-117fd7e60a32
md" #### Code 6.26"

# ╔═╡ 0ac47865-b555-4806-bab6-ede18eb823bc
begin
	Random.seed!(6)
	U = 2 .* rand(Bernoulli(), N) .- 1
	G = rand(Normal(), N)
	P = rand(MvNormal(@. b_GP*G + b_U*U))
	C = rand(MvNormal(@. b_PC*P + b_GC*G + b_U*U))
	d = DataFrame(:C => C, :P => P, :G => G, :U => U);
end

# ╔═╡ 8702f134-5470-4572-9623-a84d4dd2634d
md" #### Code 6.27"

# ╔═╡ b78b9eb4-bab3-4991-be33-26d9b7085ce3
@model function model_m6_11(P, G, C)
    a ~ Normal()
    b_PC ~ Normal()
    b_GC ~ Normal()
    μ = @. a + b_PC*P + b_GC*G
    σ ~ Exponential(1)
    C ~ MvNormal(μ, σ)
end

# ╔═╡ fe82b5ed-4836-4b83-9a14-c76deb437532
begin
	m6_11 = sample(model_m6_11(d.P, d.G, d.C), NUTS(), 1000)
	m6_11_df = DataFrame(m6_11)
	PRECIS(m6_11_df)
end

# ╔═╡ f2b666f1-dc22-465b-80ad-5d3e80e1b880
md" #### Code 6.28"

# ╔═╡ 81b3190d-c7c6-4c27-978d-0b3091dd35d4
@model function model_m6_12(P, G, U, C)
    a ~ Normal()
    b_PC ~ Normal()
    b_GC ~ Normal()
    b_U ~ Normal()
    μ = @. a + b_PC*P + b_GC*G + b_U*U
    σ ~ Exponential(1)
    C ~ MvNormal(μ, σ)
end

# ╔═╡ b798cfca-3286-4242-b0a7-bdeea06d8835
begin
	m6_12 = sample(model_m6_12(d.P, d.G, d.U, d.C), NUTS(), 1000)
	m6_12_df = DataFrame(m6_12)
	PRECIS(m6_12_df)
end

# ╔═╡ 0a2b01b3-e846-401a-8cb5-95f31d5c296a
md" ## 6.4 Confronting confounding"

# ╔═╡ 58a9e8c1-c857-4e9c-9667-6b4b45a192ab
# Code 6.29

dag_61 = Dagitty.DAG(
    :X => :Y,
    :U => :X, :A => :U,
    :A => :C, :C => :Y,
    :U => :B, :C => :B,
)

# ╔═╡ f834cb25-1de8-4a5e-8c0c-b5fc247d5592
Dagitty.all_backdoor_adjustment_sets(dag_61, :X, :Y)

# ╔═╡ be2b6b8c-8c35-469d-a4ce-bd699358302d
# Code 6.30

dag_62 = Dagitty.DAG(
    :A => :D,
    :A => :M, :M => :D,
    :S => :A, :S => :M,
    :S => :W, :W => :D,
)

# ╔═╡ 93c85c80-cfef-4e69-bceb-ab0e61398e77
all_backdoor_adjustment_sets(dag_62, :W, :D)

# ╔═╡ 83144746-ac55-48c8-aa46-e060103e86e8
implied_conditional_independencies_min(dag_62)

# ╔═╡ Cell order:
# ╟─579065ea-bc4c-49a3-a3d9-673e49071dfe
# ╠═648d7e4e-5575-4867-814f-2972fe6c6729
# ╠═3ae8474f-90d0-4653-a336-0803f7705462
# ╠═8da6e45a-9d78-479a-a8bb-57641e9a7375
# ╟─299eb270-63f8-496f-a235-400f3859babd
# ╠═40d57102-ed75-4da4-bdb1-cc18243ea570
# ╠═399bb4a5-23f2-4c9e-bae0-61ae8c668bed
# ╟─008c3130-6d49-45dd-a1da-c125c624b23b
# ╟─54c70b49-beee-4ac3-92f0-d7cd2816a6be
# ╠═0ef7a233-2a4e-4bc1-bd46-7387165b5be7
# ╟─63369b75-a94c-4365-8171-79c72051eaa1
# ╠═318c8900-a52d-4338-98d0-9dc2f1073cea
# ╠═40e28711-0a2e-40d2-a584-ae22f4de621a
# ╟─999c880d-596e-4ad2-ab4e-f2cc364ac99c
# ╠═d9bd103f-2cc5-48cc-89e0-99c08d2897bf
# ╟─aec96c68-fa64-4548-b68f-b6b0f2406182
# ╠═f4eecae9-0374-4845-b673-c11bdedc2819
# ╟─c7988aad-113e-4df4-95ed-90acfbaf60a2
# ╠═bda99113-bd50-499f-aa7c-e6eea69d48da
# ╟─ea3ecd44-ccfd-41d8-91b9-10301a49523f
# ╠═b02ec264-69ad-4d2c-a847-a9db6396631c
# ╠═04be25da-56cd-49c9-ae84-466f6fe672d8
# ╠═5a81fd3c-fc60-47a5-9da9-2e40eeb30251
# ╟─0e9cd9a9-89b3-465a-8936-f10ef02eb4d0
# ╠═4e8ea51e-dedb-4f21-8661-4220b6a6227b
# ╟─abfee492-d400-4543-a7a4-a3382bff4087
# ╠═a9b652e5-3376-43c0-b962-d905810c82e8
# ╠═020f1bbd-3f89-4256-90c0-b595712128bc
# ╠═d17704eb-1212-42a3-bcc4-f59bf3231333
# ╠═badb90cf-aaa4-4e69-92da-54ddf020b473
# ╠═1755c963-e12b-4f5a-b670-51e7cf3257f9
# ╠═fcc55095-6ae2-40eb-81bd-639a41550f0a
# ╟─61e1e17d-1347-44d6-95b9-02ecd9f58450
# ╠═b6208f9a-b80a-4b0a-8833-7a3154b0b96c
# ╠═cf869b8f-2c35-4ba3-9546-2d15cbdaf6e0
# ╠═8e1f339c-cf39-4b73-b0aa-fbd37c73e81d
# ╠═932dc19e-080f-451e-88a7-653f2ba82b1d
# ╟─a2b8eeb5-67df-4023-8631-6b001feb20c7
# ╠═7d79c5ed-e6b6-4635-9cfa-8847bd2da8ae
# ╠═fc0f6c2a-e625-4fa6-bcf7-1bc64f617c27
# ╟─2dc14c15-86b2-4b78-965a-e63fac38db6f
# ╟─03e1487b-321f-41a8-b2e8-a27d53b29fa6
# ╠═b89b0a75-3dc3-443e-af41-9b17df5dd649
# ╟─83264c24-cdbf-44a7-b194-6f4380ce3a12
# ╠═2e9c86d6-33af-4fed-ba76-304b5e8ba6fe
# ╟─ab057169-5969-4c80-8028-61dabf9b1f1a
# ╠═4005549c-b251-4485-9107-c5bfee420bdf
# ╠═ef69951e-cfdc-487b-94d6-4b21217b4a8c
# ╟─52f0e6aa-f55e-493c-b8a1-f681170b5ef8
# ╠═f8b2e102-df8a-4198-a242-58ab278d6ccc
# ╠═43526e84-61d1-40ae-8146-416b9406e7e4
# ╟─fc43c304-9e82-4e4f-b56a-44cd7f664520
# ╠═b507451b-b0e9-4786-80ab-fdae9c86f039
# ╠═7faa03d8-23b4-4792-9b08-080f357be27a
# ╟─08f49590-0dd0-40d4-ae07-74b90d923811
# ╠═052f4881-2db3-49d9-817f-c09c7abd7b80
# ╟─13242e24-b521-4181-9d57-d9dabe94c322
# ╠═b79ebc48-f738-491e-8fe2-822c1828ad90
# ╟─8902c6c3-f99a-48b5-bbd9-5ce65dda504f
# ╠═695028d4-42a9-469b-94ae-0456268e4c62
# ╠═0ce498fd-cec5-4028-872c-88f8941a0e87
# ╠═fe3d9e94-1ac1-48cf-9a58-ec409324f7f4
# ╠═3659ca0e-5ce4-4227-a6e1-b4f250e6ba30
# ╠═be9f02e7-3704-4764-bc1f-0adef06f2836
# ╠═d0c867a1-e69e-4034-9dbb-80f87dc33a2a
# ╠═a0194048-1517-4dcf-84af-dafbc0a0222c
# ╟─e1269e2b-0f88-4137-b182-2f0c8f79f6b3
# ╠═cbef9b36-bf1d-40c9-83d0-72879c683010
# ╟─0ad131a0-212a-4169-9a69-033d02cba337
# ╠═576750fa-4ff7-477b-8a17-2707c0b9b4c1
# ╠═a1832643-71b9-4728-ad53-cca2fae2da07
# ╠═b8dadb5d-bc7a-4f35-8f7e-8eb803cc6c46
# ╟─f676cc8f-8ec3-4e30-9a6d-493b65b894be
# ╠═791d6baf-d0bb-4c8f-997c-60e9d5e73e51
# ╠═a33bd041-06d5-4cd4-a88c-38bc8ff4c0cc
# ╟─8836c691-492a-4f76-9193-05389b69ba34
# ╠═c81ca3a3-8ea5-4ccb-a420-75bb0970d5b3
# ╟─b43cbceb-79be-4f66-9a4b-117fd7e60a32
# ╠═0ac47865-b555-4806-bab6-ede18eb823bc
# ╟─8702f134-5470-4572-9623-a84d4dd2634d
# ╠═b78b9eb4-bab3-4991-be33-26d9b7085ce3
# ╠═fe82b5ed-4836-4b83-9a14-c76deb437532
# ╠═f2b666f1-dc22-465b-80ad-5d3e80e1b880
# ╠═81b3190d-c7c6-4c27-978d-0b3091dd35d4
# ╠═b798cfca-3286-4242-b0a7-bdeea06d8835
# ╟─0a2b01b3-e846-401a-8cb5-95f31d5c296a
# ╠═58a9e8c1-c857-4e9c-9667-6b4b45a192ab
# ╠═f834cb25-1de8-4a5e-8c0c-b5fc247d5592
# ╠═be2b6b8c-8c35-469d-a4ce-bd699358302d
# ╠═93c85c80-cfef-4e69-bceb-ab0e61398e77
# ╠═83144746-ac55-48c8-aa46-e060103e86e8
