### A Pluto.jl notebook ###
# v0.19.3

using Markdown
using InteractiveUtils

# ╔═╡ e7613140-2875-45d0-a7b9-05019e0e55e9
using Pkg, DrWatson

# ╔═╡ ffd99df2-d429-4355-ae71-cfc2906f72d3
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
	using Plots.PlotMeasures
	using StatsBase
	using FreqTables
	using Logging
end

# ╔═╡ 57ed207c-b572-4688-b2ef-1e3c574fcea0
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

# ╔═╡ 45d98b5a-d248-41cd-9c5c-e792444dd7d9
begin
	default(label=false);
	Logging.disable_logging(Logging.Warn);
end;

# ╔═╡ 8147ea5b-1853-49bb-952d-231bc81928c5
md" ## 13.1 Example: multilevel tadpoles"

# ╔═╡ 9941ca6b-98be-4ccd-848d-6bc9e2cfd4a2
md" #### Code 13.1"

# ╔═╡ b1733234-0b5a-4868-913f-082737c16584
begin
	frogs = CSV.read(sr_datadir("reedfrogs.csv"), DataFrame)
	frogs.tank = 1:nrow(frogs)
	describe(frogs)
end

# ╔═╡ 58b836dd-cf22-4d80-a3db-94f2b405874d
md" #### Code 13.2"

# ╔═╡ a4f381c5-cf17-481f-8f42-eb9e41df61a7
@model function m13_1(S, N, tank)
    tank_size = length(levels(tank))
    a ~ filldist(Normal(0, 1.5), tank_size)
    p = logistic.(a)
    @. S ~ Binomial(N, p)
end

# ╔═╡ be6041bc-c7f1-44a8-8ca0-66fd6fdb252a
begin
	Random.seed!(1)
	m13_1_ch = sample(m13_1(frogs.surv, frogs.density, frogs.tank), NUTS(200, 0.65, init_ϵ=0.5), 1000)
	m13_1_df = DataFrame(m13_1_ch);
end

# ╔═╡ a1144feb-1264-4429-b019-d2580f1dee8e
md" #### Code 13.3"

# ╔═╡ 791f6424-bebc-4e58-adfc-c8f1ef87d4a2
md" #### Code 13.4"

# ╔═╡ 794e96b7-6a8b-414a-b7f6-0341a9cbfbac
link_fun = (r, dr) -> begin
    a = get(r, "a[$(dr.tank)]", 0)
    p = logistic(a)
    binomlogpdf(dr.density, p, dr.surv)
end

# ╔═╡ 277bef7f-964b-4073-8cbc-e098c5541843
md" #### Code 13.5"

# ╔═╡ f75d2093-821b-4d3e-9b94-7ff77240176a
md" #### Code 13.6"

# ╔═╡ 85c7d13d-7b34-42d0-a870-be47428c5a15
md" ## 13.2 Varying effects and the underfitting/overfitting trade-off."

# ╔═╡ 676f183f-fd44-47ce-a53a-b315dbdff5e1
md" #### Code 13.7"

# ╔═╡ e33a498b-395f-4684-ba2c-a39b66b91255
begin
	ā = 1.5
	σ = 1.5
	nponds = 60
	Ni = repeat([5, 10, 25, 35], inner=15);
end;

# ╔═╡ 32bc923c-fb78-4f47-8e3a-c6efae58163c
@model function m13_2(S, N, tank)
    tank_size = length(levels(tank))
    σ ~ Exponential()
    ā ~ Normal(0, 1.5)
    a ~ filldist(Normal(ā, σ), tank_size)
    p = logistic.(a)
    @. S ~ Binomial(N, p)
end

# ╔═╡ fde815d6-beb9-4114-82a0-01b8d8aadeae
begin
	Random.seed!(1)
	m13_2_ch = sample(m13_2(frogs.surv, frogs.density, frogs.tank), NUTS(200, 0.65, init_ϵ=0.2), 1000)
	m13_2_df = DataFrame(m13_2_ch);
end

# ╔═╡ 43eb9d81-f034-4c2c-bd37-619705342497
let
	m1_ll = link(m13_1_df, link_fun, eachrow(frogs))
	m1_ll = hcat(m1_ll...);
	
	m2_ll = link(m13_2_df, link_fun, eachrow(frogs))
	m2_ll = hcat(m2_ll...);
	
	compare([m1_ll, m2_ll], :waic, mnames=["m13.1", "m13.2"])
end

# ╔═╡ ccec6027-9f9d-44d1-9e91-f5f74a188acb
let
	Random.seed!()
	post = sample(resetrange(m13_2_ch), 10000)
	global post_df = DataFrame(post)
	
	propsurv_est = [
	    logistic(mean(post_df[:,"a[$i]"]))
	    for i ∈ 1:nrow(frogs)
	]
	
	scatter(propsurv_est, mc=:white, xlab="tank", ylab="proportion survival", ylim=(-0.05, 1.05))
	scatter!(frogs.propsurv, mc=:blue, ms=3)
	hline!([mean(logistic.(post_df.ā))], ls=:dash, c=:black)
	vline!([16.5, 32.5], c=:black)
	annotate!([
	        (8, 0, ("small tanks", 10)),
	        (16+8, 0, ("medium tanks", 10)),
	        (32+8, 0, ("large tanks", 10))
	    ])
end

# ╔═╡ c96fc5cc-5b3e-4032-9c7e-b43f486242e6
let
	p1 = plot(xlim=(-3, 4), xlab="log-odds survive", ylab="Density")
	
	for r ∈ first(eachrow(post_df), 100)
	    plot!(Normal(r.ā, r.σ), c=:black, alpha=0.2)
	end
	
	sim_tanks = @. rand(Normal(post_df.ā[1:8000], post_df.σ[1:8000]));
	p2 = plot(xlab="probability survive", ylab="Density", xlim=(-0.1, 1.1))
	density!(logistic.(sim_tanks), lw=2)
	
	plot(p1, p2, size=(800, 400), margin=2mm)
end

# ╔═╡ 1b684dbb-17f9-441f-b54e-371d0bf4af5e
md" #### Code 13.8-13.10"

# ╔═╡ 6ecb0714-28f3-42ac-a2c3-8a3506f801e6
let
	Random.seed!(5005)
	a_pond = rand(Normal(ā, σ), nponds);
	global dsim = DataFrame(pond=1:nponds, Ni=Ni, true_a=a_pond);
	
	# Doesn't make much sense in Julia, but anyways
	
	typeof(1:3), typeof([1,2,3])
	
	Random.seed!(1)
	dsim.Si = @. rand(Binomial(dsim.Ni, logistic(dsim.true_a)))
	dsim.p_nopool = dsim.Si ./ dsim.Ni;
	PRECIS(dsim)
end

# ╔═╡ 7ef78fb9-040a-4622-90e0-a173a6721b26
md" #### Code 13.13"

# ╔═╡ d3660350-5972-4a04-9991-2afe1d7534c0
@model function m13_3(Si, Ni, pond)
    σ ~ Exponential()
    ā ~ Normal(0, 1.5)
    a_pond ~ filldist(Normal(ā, σ), nponds)
    p = logistic.(a_pond)
    @. Si ~ Binomial(Ni, p)
end

# ╔═╡ 9d1dd75d-68fa-41e8-9a8b-5014a27061a6
begin
	Random.seed!(1)
	m13_3_ch = sample(m13_3(dsim.Si, dsim.Ni, dsim.pond), NUTS(), 1000)
	m13_3_df = DataFrame(m13_3_ch);
end

# ╔═╡ beb669bc-d900-4785-89ab-0b1f5f9def46
md" #### Code 13.14"

# ╔═╡ 338e2424-1f47-4232-97d5-f80c72c8f162
PRECIS(m13_3_df)

# ╔═╡ 94cdfc4a-c799-4668-928a-0aa265fb4f53
md" #### Code 13.15"

# ╔═╡ 7c662cb8-ee2f-4dc2-980f-e538635008df
dsim.p_partpool = [
    mean(logistic.(m13_3_df[:,"a_pond[$i]"]))
    for i ∈ 1:nponds
];

# ╔═╡ 91e82d6a-2f2a-4399-97f6-523e17d9ab77
md" #### Code 13.16"

# ╔═╡ bb520798-7ab9-4e52-8f33-aef9219bcae4
dsim.p_true = logistic.(dsim.true_a);

# ╔═╡ 79e07adc-3378-484b-8b65-45275741ca5b
md" #### Code 13.17 - 13.19"

# ╔═╡ d28addae-fc36-4822-9ffc-bb831a0ad03d
begin
	nopool_error = @. abs(dsim.p_nopool - dsim.p_true)
	partpool_error = @. abs(dsim.p_partpool - dsim.p_true)
end

# ╔═╡ 89fc29e9-9932-45c9-ba71-ddea749a3eb1
let
	scatter(nopool_error, xlab="pond", ylab="absolute error")
	scatter!(partpool_error, mc=:white)
end

# ╔═╡ fc227622-466f-4e6f-8b10-c9062915704d
md" #### Code 13.19 - 13.20"

# ╔═╡ 31a10846-a373-4e38-87e8-7bed82ae8707
let
	dsim.nopool_error = nopool_error
	dsim.partpool_error = partpool_error
	
	gb = groupby(dsim, :Ni)
	nopool_avg = combine(gb, :nopool_error => mean)
	partpool_avg = combine(gb, :partpool_error => mean);
	nopool_avg, partpool_avg

	ā = 1.5
	σ = 1.5
	nponds = 60
	Ni = repeat([5, 10, 25, 35], inner=15)
	a_pond = rand(Normal(ā, σ), nponds)
	
	dsim = DataFrame(pond=1:nponds, Ni=Ni, true_a=a_pond)
	dsim.Si = @. rand(Binomial(dsim.Ni, logistic(dsim.true_a)))
	dsim.p_nopool = dsim.Si ./ dsim.Ni
	
	m13_3_ch = sample(m13_3(dsim.Si, dsim.Ni, dsim.pond), NUTS(), 1000)
	m13_3_df = DataFrame(m13_3_ch)

	dsim.p_partpool = [
	    mean(logistic.(m13_3_df[:,"a_pond[$i]"]))
	    for i ∈ 1:nponds
	]
	dsim.p_true = logistic.(dsim.true_a)
	nopool_error = @. abs(dsim.p_nopool - dsim.p_true)
	partpool_error = @. abs(dsim.p_partpool - dsim.p_true)
	
	scatter(nopool_error, xlab="pond", ylab="absolute error")
	scatter!(partpool_error, mc=:white)
end

# ╔═╡ 1f195800-e327-46cf-8bf5-f728756f0132
md" ## 13.3 More than one type of cluster."

# ╔═╡ d66a906b-d034-4144-bd4f-d1c3dff7015d
md" #### Code 13.21"

# ╔═╡ 9293a32f-4296-4f96-88ac-00ac51c23360
begin
	chimpanzees = CSV.read(sr_datadir("chimpanzees.csv"), DataFrame)
	chimpanzees.treatment = 1 .+ chimpanzees.prosoc_left .+ 2*chimpanzees.condition;
end;

# ╔═╡ 03f31f7c-b6fa-46b3-9785-2b70abbe9dd1
@model function m13_4(pulled_left, actor, block_id, treatment)
    σ_a ~ Exponential()
    σ_g ~ Exponential()
    ā ~ Normal(0, 1.5)
    actors_count = length(levels(actor))
    blocks_count = length(levels(block_id))
    treats_count = length(levels(treatment))
    a ~ filldist(Normal(ā, σ_a), actors_count)
    g ~ filldist(Normal(0, σ_g), blocks_count)
    b ~ filldist(Normal(0, 0.5), treats_count)
    
    p = @. logistic(a[actor] + g[block_id] + b[treatment])
    @. pulled_left ~ Binomial(1, p)
end

# ╔═╡ 012f1b59-ddb2-4906-b5fb-cad11527af97
begin
	Random.seed!(13)
	m13_4_ch = sample(m13_4(chimpanzees.pulled_left, chimpanzees.actor, chimpanzees.block, chimpanzees.treatment),
		NUTS(), 4000)
	m13_4_df = DataFrame(m13_4_ch);
end

# ╔═╡ 448d4c51-ae2a-4154-9199-756114f6eecf
md" #### Code 13.22"

# ╔═╡ c061efe3-bbfe-4736-8f97-091ec5088efc
PRECIS(m13_4_df)

# ╔═╡ be23b9da-4de5-4388-9a64-fedbccaae986
coeftab_plot(m13_4_df, size=(800,600))

# ╔═╡ 678a47cb-6d8b-4994-a713-cc8c274f7411
md" #### Code 13.23"

# ╔═╡ a621a173-458f-4872-93f7-4e57e19856b9
@model function m13_5(pulled_left, actor, treatment)
    σ_a ~ Exponential()
    ā ~ Normal(0, 1.5)
    actors_count = length(levels(actor))
    treats_count = length(levels(treatment))
    a ~ filldist(Normal(ā, σ_a), actors_count)
    b ~ filldist(Normal(0, 0.5), treats_count)
    
    p = @. logistic(a[actor] + b[treatment])
    @. pulled_left ~ Binomial(1, p)
end

# ╔═╡ 9bef0f5a-28e5-4074-b9d9-f7e52e452fb6
begin
	Random.seed!(14)
	m13_5_ch = sample(m13_5(chimpanzees.pulled_left, chimpanzees.actor, chimpanzees.treatment), NUTS(), 4000)
	m13_5_df = DataFrame(m13_5_ch);
end

# ╔═╡ 65beeb3d-d0f4-4931-ac20-b7b22be19eb2
md" #### Code 13.24"

# ╔═╡ 44783098-f38e-44fe-86e1-946d02f3f9db
l_fun1 = (r, dr) -> begin
    a = get(r, "a[$(dr.actor)]", 0)
    g = get(r, "g[$(dr.block)]", 0)
    b = get(r, "b[$(dr.treatment)]", 0)
    p = logistic(a + g + b)
    binomlogpdf(1, p, dr.pulled_left)
end

# ╔═╡ 73c08bc2-0ebe-4f5a-a769-72d41539a6a0
begin
	m13_4_ll = link(m13_4_df, l_fun1, eachrow(chimpanzees))
	m13_4_ll = hcat(m13_4_ll...)
end

# ╔═╡ 0e7f2b1b-4ddd-44e9-9c99-4c9e604eb89a
l_fun2 = (r, dr) -> begin
    a = get(r, "a[$(dr.actor)]", 0)
    b = get(r, "b[$(dr.treatment)]", 0)
    p = logistic(a + b)
    binomlogpdf(1, p, dr.pulled_left)
end

# ╔═╡ 3c5564d4-7be2-4028-bb15-ba0313da6782
begin
	m13_5_ll = link(m13_5_df, l_fun2, eachrow(chimpanzees))
	m13_5_ll = hcat(m13_5_ll...);
end

# ╔═╡ 51067681-8aa3-471d-85df-e78786cec73d
compare([m13_4_ll, m13_5_ll], :waic, mnames=["m4", "m5"])

# ╔═╡ 55780648-3d5e-449b-94cd-312b5d70dad6
md" #### Code 13.25"

# ╔═╡ 73e31e16-bbe9-4c9a-96bb-194c8110deba
@model function m13_6(pulled_left, actor, block_id, treatment)
    σ_a ~ Exponential()
    σ_g ~ Exponential()
    σ_b ~ Exponential()
    ā ~ Normal(0, 1.5)
    actors_count = length(levels(actor))
    blocks_count = length(levels(block_id))
    treats_count = length(levels(treatment))
    a ~ filldist(Normal(ā, σ_a), actors_count)
    g ~ filldist(Normal(0, σ_g), blocks_count)
    b ~ filldist(Normal(0, σ_b), treats_count)
    
    p = @. logistic(a[actor] + g[block_id] + b[treatment])
    @. pulled_left ~ Binomial(1, p)
end

# ╔═╡ b2598a77-4d7b-4878-b555-d1afa47e0eb6
begin
	Random.seed!(15)
	m13_6_ch = sample(m13_6(chimpanzees.pulled_left, chimpanzees.actor, chimpanzees.block, chimpanzees.treatment), NUTS(), 4000)
	m13_6_df = DataFrame(m13_6_ch);
end

# ╔═╡ e5df33fd-885f-419f-86ca-de7233b901e4
PRECIS(m13_4_df[:,r"b"])

# ╔═╡ 95164e3b-60db-4e7f-84ed-61d1b2978ff0
PRECIS(m13_6_df[:,r"b"])

# ╔═╡ a494c2d9-e049-4fbe-aaa6-202810673ace
md" ## 13.4 Divergent transitions and non-centered priors."

# ╔═╡ 0eecba9e-0563-4eba-9b2f-679cf1140fb8
md" #### Code 13.26"

# ╔═╡ b1105cc8-bbbc-4aea-ad2f-c12f238e55ac
@model function m13_7(N)
    v ~ Normal(0, 3)
    x ~ Normal(0, exp(v))
end

# ╔═╡ 1acd31f3-b759-454e-aafe-9b36073f5f4f
begin
	Random.seed!(5)
	m13_7_ch = sample(m13_7(1), NUTS(), 1000)
end

# ╔═╡ f6c7c0dd-eba5-4118-b854-a175ddde0ad5
md" #### Code 13.27"

# ╔═╡ fe449e67-fe7b-4713-a84b-f4a91f0d9380
@model function m13_7nc(N)
    v ~ Normal(0, 3)
    z ~ Normal()
    x = z * exp(v)
end

# ╔═╡ 8b56d961-3cca-4431-8d97-c33bd761387c
begin
	Random.seed!(5)
	m13_7nc_ch = sample(m13_7nc(1), NUTS(), 1000)
end

# ╔═╡ a52b6a55-cf23-465b-b7d9-0bb741372d40
md" #### Code 13.28"

# ╔═╡ d89aecb5-5fa1-4d94-b216-52a101c87747
md"

!!! note

There is no way to get amount of divergent samples, but they could be estimated by comparing `ess` values from the chain."

# ╔═╡ c38dc9ee-5cbb-47dd-a14a-cd185b4b4c46
begin
	Random.seed!(13)
	m13_4b_ch = sample(m13_4(chimpanzees.pulled_left, chimpanzees.actor, chimpanzees.block, chimpanzees.treatment),
		NUTS(0.95, init_ϵ=0.1), 4000)
end

# ╔═╡ f05cec5b-bb78-4e26-b5ea-e3e51c42a793
let
	t = ess_rhat(m13_4_ch)
	ess_4 = t[:,:ess]
	t = ess_rhat(m13_4b_ch)
	ess_4b = t[:,:ess]
	
	plot(ess_4, lw=2, label="ESS m13_4")
	plot!(ess_4b, lw=2, label="ESS m13_4b")
end

# ╔═╡ 9b59ade6-affb-4a20-aacb-c02967c201b0
md" #### Code 13.29"

# ╔═╡ 486e7146-4bd2-4465-a562-bd2052c123e7
@model function m13_4nc(pulled_left, actor, block_id, treatment)
    σ_a ~ Exponential()
    σ_g ~ Exponential()
    ā ~ Normal(0, 1.5)
    actors_count = length(levels(actor))
    blocks_count = length(levels(block_id))
    treats_count = length(levels(treatment))
    b ~ filldist(Normal(0, 0.5), treats_count)
    z ~ filldist(Normal(), actors_count)
    x ~ filldist(Normal(), blocks_count)
    a = @. ā + σ_a*z
    g = σ_g*x
    
    p = @. logistic(a[actor] + g[block_id] + b[treatment])
    @. pulled_left ~ Binomial(1, p)
end

# ╔═╡ 12e79119-c060-4ac4-92e2-9658a1e9d587
begin
	Random.seed!(13)
	m13_4nc_ch = sample(m13_4nc(chimpanzees.pulled_left, chimpanzees.actor, chimpanzees.block, chimpanzees.treatment),
		NUTS(), 4000);
end

# ╔═╡ d6b6a49f-a4d7-473b-9e46-324cd335b80a
md" #### Code 13.30"

# ╔═╡ a4995532-cb86-4c30-bd3d-46b8f6a9253f
let
	t = ess_rhat(m13_4_ch)
	ess_4 = t[:,:ess]
	t = ess_rhat(m13_4nc_ch)
	ess_4nc = t[:,:ess]
	
	lims = extrema(vcat(ess_4, ess_4nc)) .+ (-100, 100)
	plot(xlim=lims, ylims=lims, xlab="n_eff (centered)", ylab="n_eff (non-centered)", size=(500,500))
	scatter!(ess_4, ess_4nc)
	plot!(identity, c=:gray, s=:dash)
end

# ╔═╡ ef1c6253-944e-4647-aed2-c71df9ec2906
md" ## 13.5 Multilevel posterior predictions."

# ╔═╡ 9bec6933-dd24-40b5-92a7-ca8956616ddd
md" #### Code 13.31"

# ╔═╡ 002ac714-3e07-40fb-83a2-a2f1c59cfc0f
let
	chimp = 2
	d_pred = DataFrame(
	    actor = fill(chimp, 4),
	    treatment = 1:4,
	    block = fill(1, 4)
	)
	
	l_fun = (r, dr) -> begin
	    a = get(r, "a[$(dr.actor)]", 0)
	    g = get(r, "g[$(dr.block)]", 0)
	    b = get(r, "b[$(dr.treatment)]", 0)
	    logistic(a + g + b)
	end
	
	p = link(m13_4_df, l_fun, eachrow(d_pred))
	p = hcat(p...)
	p_μ = mean.(eachcol(p))
	p_ci = PI.(eachcol(p));
end

# ╔═╡ 9384e567-cc37-425e-928f-9ec8cd8ec8e5
md" #### Code 13.32"

# ╔═╡ 14745152-0313-489e-a668-b17cfa7d0cfd
begin
	post13_4 = sample(resetrange(m13_4_ch), 2000)
	post13_4_df = DataFrame(post13_4)
	describe(post13_4_df)
end

# ╔═╡ fe42ac20-c7f2-416b-af21-e4dccbc4f436
md" #### Code 13.33"

# ╔═╡ aad99eb9-7ee2-4beb-b461-ccdc45ba04b9
density(post13_4_df."a[5]")

# ╔═╡ 16edfd3e-26e0-47b7-a36c-e4144ca53d38
md" #### Code 13.34"

# ╔═╡ 30b9de2d-97e5-4d7c-9042-1018dfc53b3a
p_link = (treatment, actor, block_id) -> begin
    logodds = 
        getproperty(post13_4_df, "a[$actor]") + 
        getproperty(post13_4_df, "b[$block_id]") + 
        getproperty(post13_4_df, "g[$treatment]")
    logistic.(logodds)
end

# ╔═╡ 43f08753-8307-4094-b183-5ef7017e3a1e
md" #### Code 13.35"

# ╔═╡ e90c9e4c-5c00-4253-b6b8-ad9c7b749e95
begin
	p_raw = p_link.(1:4, 2, 1)
	p_raw = hcat(p_raw...)
	p_μ = mean.(eachcol(p_raw))
	p_ci = PI.(eachcol(p_raw));
end

# ╔═╡ b34a1f96-0d66-4031-9a99-20e73780c7c5
md" #### Code 13.36"

# ╔═╡ c57afb19-d62d-4e8c-b2dd-2b176a6280e0
p_link_abar = treatment -> begin
    logodds = post13_4_df.ā + getproperty(post13_4_df, "b[$treatment]")
    logistic.(logodds)
end

# ╔═╡ d015e9cf-a420-43eb-a058-b5eac0e10662
md" #### Code 13.37"

# ╔═╡ d359891c-eefd-4e2a-94c7-0d20f35c61f6
let
	p_raw = p_link_abar.(1:4)
	p_raw = hcat(p_raw...)
	p_μ = mean.(eachcol(p_raw))
	p_ci = PI.(eachcol(p_raw))
	p_ci = vcat(p_ci'...)
	
	plot(xlab="treatment", ylab="proportion pulled left", title="average actor", ylim=(0, 1))
	plot!(["R/N", "L/N", "R/P", "L/P"], [p_μ p_μ], fillrange=p_ci, fillalpha=0.2, c=:black, lw=1.5)
end

# ╔═╡ d4de4418-316d-4976-a5aa-0923346eb0be
md" #### Code 13.38"

# ╔═╡ 7f4622a8-f2ea-41bf-a04e-f9f69faec4a3
let
	Random.seed!(1)
	a_sim = rand.(Normal.(post13_4_df.ā, post13_4_df.σ_a))
	
	p_link_asim = treatment -> begin
	    logodds = a_sim + getproperty(post13_4_df, "b[$treatment]")
	    logistic.(logodds)
	end
	
	global p_raw_asim = p_link_asim.(1:4)
	p_raw_asim = hcat(p_raw_asim...)
	p_μ = mean.(eachcol(p_raw_asim))
	p_ci = PI.(eachcol(p_raw_asim))
	p_ci = vcat(p_ci'...)
	
	plot(xlab="treatment", ylab="proportion pulled left", title="marginal of actor", ylim=(0, 1))
	plot!(["R/N", "L/N", "R/P", "L/P"], [p_μ p_μ], fillrange=p_ci, fillalpha=0.2, c=:black, lw=1.5)
end

# ╔═╡ 41920919-0df9-42bf-9749-268a3a5906c2
md" #### Code 13.39"

# ╔═╡ 7864f790-f9a0-447b-aad5-c56c0f3e1f34
let
	p = plot(xlab="treatment", ylab="proportion pulled left", title="simulated actors", ylim=(0, 1))
	
	for r in first(eachrow(p_raw_asim), 100)
	    plot!(["R/N", "L/N", "R/P", "L/P"], r, c=:black, alpha=0.2)
	end
	p
end

# ╔═╡ Cell order:
# ╠═57ed207c-b572-4688-b2ef-1e3c574fcea0
# ╠═e7613140-2875-45d0-a7b9-05019e0e55e9
# ╠═ffd99df2-d429-4355-ae71-cfc2906f72d3
# ╠═45d98b5a-d248-41cd-9c5c-e792444dd7d9
# ╟─8147ea5b-1853-49bb-952d-231bc81928c5
# ╟─9941ca6b-98be-4ccd-848d-6bc9e2cfd4a2
# ╠═b1733234-0b5a-4868-913f-082737c16584
# ╟─58b836dd-cf22-4d80-a3db-94f2b405874d
# ╠═a4f381c5-cf17-481f-8f42-eb9e41df61a7
# ╠═be6041bc-c7f1-44a8-8ca0-66fd6fdb252a
# ╟─a1144feb-1264-4429-b019-d2580f1dee8e
# ╠═32bc923c-fb78-4f47-8e3a-c6efae58163c
# ╠═fde815d6-beb9-4114-82a0-01b8d8aadeae
# ╠═791f6424-bebc-4e58-adfc-c8f1ef87d4a2
# ╠═794e96b7-6a8b-414a-b7f6-0341a9cbfbac
# ╠═43eb9d81-f034-4c2c-bd37-619705342497
# ╟─277bef7f-964b-4073-8cbc-e098c5541843
# ╠═ccec6027-9f9d-44d1-9e91-f5f74a188acb
# ╟─f75d2093-821b-4d3e-9b94-7ff77240176a
# ╠═c96fc5cc-5b3e-4032-9c7e-b43f486242e6
# ╟─85c7d13d-7b34-42d0-a870-be47428c5a15
# ╟─676f183f-fd44-47ce-a53a-b315dbdff5e1
# ╠═e33a498b-395f-4684-ba2c-a39b66b91255
# ╠═1b684dbb-17f9-441f-b54e-371d0bf4af5e
# ╠═6ecb0714-28f3-42ac-a2c3-8a3506f801e6
# ╟─7ef78fb9-040a-4622-90e0-a173a6721b26
# ╠═d3660350-5972-4a04-9991-2afe1d7534c0
# ╠═9d1dd75d-68fa-41e8-9a8b-5014a27061a6
# ╠═beb669bc-d900-4785-89ab-0b1f5f9def46
# ╠═338e2424-1f47-4232-97d5-f80c72c8f162
# ╟─94cdfc4a-c799-4668-928a-0aa265fb4f53
# ╠═7c662cb8-ee2f-4dc2-980f-e538635008df
# ╟─91e82d6a-2f2a-4399-97f6-523e17d9ab77
# ╠═bb520798-7ab9-4e52-8f33-aef9219bcae4
# ╟─79e07adc-3378-484b-8b65-45275741ca5b
# ╠═d28addae-fc36-4822-9ffc-bb831a0ad03d
# ╠═89fc29e9-9932-45c9-ba71-ddea749a3eb1
# ╠═fc227622-466f-4e6f-8b10-c9062915704d
# ╠═31a10846-a373-4e38-87e8-7bed82ae8707
# ╟─1f195800-e327-46cf-8bf5-f728756f0132
# ╟─d66a906b-d034-4144-bd4f-d1c3dff7015d
# ╠═9293a32f-4296-4f96-88ac-00ac51c23360
# ╠═03f31f7c-b6fa-46b3-9785-2b70abbe9dd1
# ╠═012f1b59-ddb2-4906-b5fb-cad11527af97
# ╟─448d4c51-ae2a-4154-9199-756114f6eecf
# ╠═c061efe3-bbfe-4736-8f97-091ec5088efc
# ╠═be23b9da-4de5-4388-9a64-fedbccaae986
# ╟─678a47cb-6d8b-4994-a713-cc8c274f7411
# ╠═a621a173-458f-4872-93f7-4e57e19856b9
# ╠═9bef0f5a-28e5-4074-b9d9-f7e52e452fb6
# ╟─65beeb3d-d0f4-4931-ac20-b7b22be19eb2
# ╠═44783098-f38e-44fe-86e1-946d02f3f9db
# ╠═73c08bc2-0ebe-4f5a-a769-72d41539a6a0
# ╠═0e7f2b1b-4ddd-44e9-9c99-4c9e604eb89a
# ╠═3c5564d4-7be2-4028-bb15-ba0313da6782
# ╠═51067681-8aa3-471d-85df-e78786cec73d
# ╠═55780648-3d5e-449b-94cd-312b5d70dad6
# ╠═73e31e16-bbe9-4c9a-96bb-194c8110deba
# ╠═b2598a77-4d7b-4878-b555-d1afa47e0eb6
# ╠═e5df33fd-885f-419f-86ca-de7233b901e4
# ╠═95164e3b-60db-4e7f-84ed-61d1b2978ff0
# ╟─a494c2d9-e049-4fbe-aaa6-202810673ace
# ╟─0eecba9e-0563-4eba-9b2f-679cf1140fb8
# ╠═b1105cc8-bbbc-4aea-ad2f-c12f238e55ac
# ╠═1acd31f3-b759-454e-aafe-9b36073f5f4f
# ╟─f6c7c0dd-eba5-4118-b854-a175ddde0ad5
# ╠═fe449e67-fe7b-4713-a84b-f4a91f0d9380
# ╠═8b56d961-3cca-4431-8d97-c33bd761387c
# ╟─a52b6a55-cf23-465b-b7d9-0bb741372d40
# ╟─d89aecb5-5fa1-4d94-b216-52a101c87747
# ╠═c38dc9ee-5cbb-47dd-a14a-cd185b4b4c46
# ╠═f05cec5b-bb78-4e26-b5ea-e3e51c42a793
# ╟─9b59ade6-affb-4a20-aacb-c02967c201b0
# ╠═486e7146-4bd2-4465-a562-bd2052c123e7
# ╠═12e79119-c060-4ac4-92e2-9658a1e9d587
# ╟─d6b6a49f-a4d7-473b-9e46-324cd335b80a
# ╠═a4995532-cb86-4c30-bd3d-46b8f6a9253f
# ╟─ef1c6253-944e-4647-aed2-c71df9ec2906
# ╟─9bec6933-dd24-40b5-92a7-ca8956616ddd
# ╠═002ac714-3e07-40fb-83a2-a2f1c59cfc0f
# ╠═9384e567-cc37-425e-928f-9ec8cd8ec8e5
# ╠═14745152-0313-489e-a668-b17cfa7d0cfd
# ╟─fe42ac20-c7f2-416b-af21-e4dccbc4f436
# ╠═aad99eb9-7ee2-4beb-b461-ccdc45ba04b9
# ╟─16edfd3e-26e0-47b7-a36c-e4144ca53d38
# ╠═30b9de2d-97e5-4d7c-9042-1018dfc53b3a
# ╟─43f08753-8307-4094-b183-5ef7017e3a1e
# ╠═e90c9e4c-5c00-4253-b6b8-ad9c7b749e95
# ╟─b34a1f96-0d66-4031-9a99-20e73780c7c5
# ╠═c57afb19-d62d-4e8c-b2dd-2b176a6280e0
# ╟─d015e9cf-a420-43eb-a058-b5eac0e10662
# ╠═d359891c-eefd-4e2a-94c7-0d20f35c61f6
# ╠═d4de4418-316d-4976-a5aa-0923346eb0be
# ╠═7f4622a8-f2ea-41bf-a04e-f9f69faec4a3
# ╟─41920919-0df9-42bf-9749-268a3a5906c2
# ╠═7864f790-f9a0-447b-aad5-c56c0f3e1f34
