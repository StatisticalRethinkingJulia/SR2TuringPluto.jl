### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 14088e14-1fe6-41bd-a70a-38595fde2d41
using Pkg, DrWatson

# ╔═╡ 10cd1a4b-7190-4a1e-a867-80555d4bf730
begin
	using Optim
	using Turing
	using DataFrames
	using CSV
	using Random
	using StatisticalRethinking
	using StatisticalRethinking: link
	using StatisticalRethinkingPlots
	using ParetoSmooth
	using StatsPlots
	using StatsBase
	using FreqTables
	using Logging
end

# ╔═╡ 6975f6c3-8eec-47d9-ade5-88e8447b2dce
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


# ╔═╡ 22f28b65-a270-4778-8933-57459e1249c7
begin
	default(label=false)
	Logging.disable_logging(Logging.Warn);
end;

# ╔═╡ f3a514c8-bc3e-487f-af11-b43ff41c34e6
md" ## 11.1 Binomial regression."

# ╔═╡ e05e454f-c5b8-4f09-b745-a0265b5718a8
md" #### Code 11.1"

# ╔═╡ a67ed97c-59ed-461b-803b-f2cb42498dda
chimpanzees = CSV.read(sr_datadir("chimpanzees.csv"), DataFrame; delim=';');

# ╔═╡ d6404e85-ceea-47c1-bc0e-c7c3deeb0ad6
md" #### Code 11.2"

# ╔═╡ dcfef783-e51f-400b-9702-9068b8ef45b7
chimpanzees[!,:treatment] = 1 .+ chimpanzees.prosoc_left .+ 2*chimpanzees.condition;

# ╔═╡ 62cf8524-0b46-434f-b813-f7befa9eb7f9
md" #### Code 11.3"

# ╔═╡ 079e0242-0d68-4e72-be44-9e317733071c
freqtable(chimpanzees, :treatment, :prosoc_left, subset=chimpanzees.condition .== 0)

# ╔═╡ e88ff706-532d-4890-ad6d-928aed3da558
freqtable(chimpanzees, :treatment, :prosoc_left, subset=chimpanzees.condition .== 1)

# ╔═╡ addadf8a-a757-4927-9e78-0b64a9bcf3eb
md" #### Code 11.4"

# ╔═╡ 247d18b1-dd98-49db-9185-2f76dae4b346
@model function ppl11_1(pulled_left)
    a ~ Normal(0, 10)
    p = logistic(a)     # inverse of the `logit` function
    pulled_left ~ Binomial(1, p)
end

# ╔═╡ c5e308ce-2a23-4773-b8d1-9dfb6a011d16
md" #### Code 11.5"

# ╔═╡ 35afca2b-8537-4ff1-a1c1-3021f3af3d4c
begin
	Random.seed!(1999)
	prior11_1 = sample(ppl11_1(chimpanzees.pulled_left), Prior(), 10000)
	prior11_1_df = DataFrame(prior11_1);
end

# ╔═╡ 38668cc3-5bcb-42fe-b89f-fed08dfacce5
md" #### Code 11.6"

# ╔═╡ 0195e6da-43f9-4ed3-bcad-7bee87fe23cf
let
	p = logistic.(prior11_1_df.a)
	density(p, bandwidth=0.01)
end

# ╔═╡ 3172be0c-e3f6-43c3-b740-137584b19bcf
md" #### Code 11.7"

# ╔═╡ 1322553d-53ba-4166-a65f-a1350158dd46
@model function ppl11_2(pulled_left, treatment)
    a ~ Normal(0, 1.5)
    treat_levels = length(levels(treatment))
    b ~ MvNormal(zeros(treat_levels), 10)
    
    p = @. logistic(a + b[treatment])
    for i ∈ eachindex(pulled_left)
        pulled_left[i] ~ Binomial(1, p[i])
    end
end

# ╔═╡ ae02befe-b6b1-49d4-8f85-00f59cd57f92
begin
	Random.seed!(1999)
	prior_chain11_2t = sample(ppl11_2(chimpanzees.pulled_left, chimpanzees.treatment), Prior(), 10000)
	prior11_2t_df = DataFrame(prior_chain11_2t)
	describe(prior11_2t_df)
end

# ╔═╡ 140245a2-38cc-440a-96f2-7e375b60ef7b
md" #### Code 11.8"

# ╔═╡ 93cdbcce-1958-4476-a32a-9b1cf3ab9f17
let
	f = i -> @. logistic(prior11_2t_df.a + prior11_2t_df[!,"b[$i]"])
	p1 = map(f, 1:4)
	density(abs.(p1[1] .- p1[2]), bandwidth=0.01)
end

# ╔═╡ b5dcad20-acc6-446b-8231-cf9d6eb6bf6d
md" #### Code 11.9"

# ╔═╡ 1430a43d-d75d-4aff-8f43-a3c5e047c6eb
@model function m11_3(pulled_left, treatment)
    a ~ Normal(0, 1.5)
    treat_levels = length(levels(treatment))
    b ~ MvNormal(zeros(treat_levels), 0.5)
    
    p = @. logistic(a + b[treatment])
    for i ∈ eachindex(pulled_left)
        pulled_left[i] ~ Binomial(1, p[i])
    end
end

# ╔═╡ a86c573b-d227-46f7-bd5b-81c44ef2b326
begin
	Random.seed!(1999)
	prior_chain3 = sample(m11_3(chimpanzees.pulled_left, chimpanzees.treatment), Prior(), 10000)
	prior3 = DataFrame(prior_chain3);
end

# ╔═╡ e172123a-8186-4038-b9d6-edd5deffefd9
let
	f = i -> @. logistic(prior3.a + prior3[!,"b[$i]"])
	p2 = map(f, 1:4)
	mean(abs.(p2[1] .- p2[2]))
end

# ╔═╡ 776b5efd-f42b-49bc-98ed-4893eb216e2c
md" #### Code 11.10"

# ╔═╡ 9780cf8c-ca3f-486e-bae3-ca480dfd8d72
dat_list = chimpanzees[!,[:pulled_left, :actor, :treatment]];

# ╔═╡ b9d8d506-d88c-48a6-9e39-d33f0984047c
md" #### Code 11.11"

# ╔═╡ ca8c5830-cd5d-421f-af9e-a0da2656d953
@model function m11_4(actor, treatment, pulled_left)
    act_levels = length(levels(actor))
    a ~ MvNormal(zeros(act_levels), 1.5)
    treat_levels = length(levels(treatment))
    b ~ MvNormal(zeros(treat_levels), 0.5)
    
    p = @. logistic(a[actor] + b[treatment])
    for i ∈ eachindex(pulled_left)
        pulled_left[i] ~ Binomial(1, p[i])
    end
end

# ╔═╡ a80e763b-9f3b-44aa-9abc-26362b4618bd
begin
	m11_4_chain = sample(m11_4(dat_list.actor, dat_list.treatment, dat_list.pulled_left), NUTS(), 1000)
	m11_4_df = DataFrame(m11_4_chain)
	describe(m11_4_df)
end

# ╔═╡ b6a5f2a8-6e45-47ba-be97-961cfeb57699
md" #### Code 11.12"

# ╔═╡ e0fb292f-8aba-4ca9-ba0c-1791cf063623
let
	p_left = DataFrame(map(i -> "V$i" => logistic.(m11_4_df[:,"a[$i]"]), 1:7)...);
	coeftab_plot(p_left)
end

# ╔═╡ c778178e-a6ab-4071-90d7-9cf543ef2606
md" #### Code 11.13"

# ╔═╡ da6198fa-c7f2-4277-98f7-d6ba2e092815
let
	names = ["R/N", "L/N", "R/P", "L/P"]
	labs = DataFrame(map(i -> names[i] => m11_4_df[:,"b[$i]"], 1:4)...)
	coeftab_plot(labs)
end

# ╔═╡ 14bc80d0-3528-4cda-9311-02af2c322fc4
md" #### Code 11.14"

# ╔═╡ 5dd96185-34b6-460e-8f7f-e34d0c5eb469
let
	diffs = DataFrame(
	    db13=m11_4_df."b[1]" .- m11_4_df."b[3]",
	    db24=m11_4_df."b[2]" .- m11_4_df."b[4]",
	)
	coeftab_plot(diffs)
end

# ╔═╡ ea6f1ad4-781a-47a6-9cee-785da25f5967
md" #### Code 11.15"

# ╔═╡ 8adfcdf2-3da5-4d30-9301-b9c2e63d5490
let
	gd = groupby(chimpanzees, [:actor, :treatment])
	c = combine(gd, :pulled_left => mean => :val)
	global pl = unstack(c, :actor, :treatment, :val)
	pl[1,:]
end

# ╔═╡ 4e297c9d-e181-41ce-b1be-d5a6669bb903
md" #### Code 11.16"

# ╔═╡ c8bd4ef6-7a94-46c7-be9a-6943043a6670
let
	names = ["R/N", "L/N", "R/P", "L/P"]
	p = plot(ylims=(0, 1.1), ylab="proportion left lever", showaxis=:y, xticks=false)
	hline!([0.5], c=:gray, s=:dash)
	for actor in 1:7
	    ofs = (actor-1)*4
	    actor > 1 && vline!([ofs+0.5], c=:gray)
	    plot!([ofs+1,ofs+3], collect(pl[actor,["1","3"]]), lw=2, m=:o, c=:blue)
	    plot!([ofs+2,ofs+4], collect(pl[actor,["2","4"]]), lw=2, m=:o, c=:blue)
	    anns = [
	        (ofs+idx, pl[actor,string(idx)]+.04, (name, 8))
	        for (idx,name) ∈ enumerate(names)
	    ]
	    actor != 2 && annotate!(anns)
	end
	
	annotate!([
	    (2.5 + (idx-1)*4, 1.1, ("actor $idx", 8))
	    for idx ∈ 1:7
	])
	p
end

# ╔═╡ e80cb9db-b93b-4b43-a55a-ee547fbc9d48
md" #### Code 11.17"

# ╔═╡ fa2b5e49-933f-4fba-886c-0cf99f81eab5
let
	l_fun = (r, (ai, bi)) -> logistic(get(r, "a[$ai]", missing) + get(r, "b[$bi]", missing))
	p_post = link(m11_4_df, l_fun, Iterators.product(1:7, 1:4))
	p_μ = map(mean, p_post)
	p_ci = map(PI, p_post);
	
	# visualize mean and cred intervals of posterior distribution
	
	# compute relative intervals 
	rel_ci = map(idx -> (p_μ[idx]-p_ci[idx][1], p_ci[idx][2]-p_μ[idx]), CartesianIndices(p_ci))
	
	names = ["R/N", "L/N", "R/P", "L/P"]
	p = plot(ylims=(0, 1.1), ylab="proportion left lever", title="Posterior", showaxis=:y, xticks=false)
	hline!([0.5], c=:gray, s=:dash)
	for actor in 1:7
	    ofs = (actor-1)*4
	    actor > 1 && vline!([ofs+0.5], c=:gray)
	    err = [rel_ci[actor,1], rel_ci[actor,3]]
	    plot!([ofs+1,ofs+3], collect(p_μ[actor,[1,3]]), err=err, lw=2, m=:o, c=:blue)
	    err = [rel_ci[actor,2], rel_ci[actor,4]]    
	    plot!([ofs+2,ofs+4], collect(p_μ[actor,[2,4]]), err=err, lw=2, m=:o, c=:blue)
	    anns = [
	        (ofs+idx, p_μ[actor,idx]+.04, (name, 8))
	        for (idx,name) ∈ enumerate(names)
	    ]
	    actor != 2 && annotate!(anns)
	end
	
	annotate!([
	    (2.5 + (idx-1)*4, 1.1, ("actor $idx", 8))
	    for idx ∈ 1:7
	])
	p
end

# ╔═╡ 47a558cd-abbb-4dcb-8130-625d433d94fa
md" #### Code 11.18"

# ╔═╡ 82744b52-b300-4821-b095-19229835c5d4
begin
	chimpanzees.side = chimpanzees.prosoc_left .+ 1
	chimpanzees.cond = chimpanzees.condition .+ 1;
end;

# ╔═╡ 9a9ca097-db59-4e94-b2d5-47f4dbde48a3
md" #### Code 11.19"

# ╔═╡ ea03ef5e-6332-41e7-a07b-55109cb97697
@model function m11_5(actor, side, cond, pulled_left)
    act_levels = length(levels(actor))
    a ~ MvNormal(zeros(act_levels), 1.5)
    side_levels = length(levels(side))
    bs ~ MvNormal(zeros(side_levels), 0.5)    
    cond_levels = length(levels(cond))
    bc ~ MvNormal(zeros(cond_levels), 0.5)
    
    p = @. logistic(a[actor] + bs[side] + bc[cond])
    for i ∈ eachindex(pulled_left)
        pulled_left[i] ~ Binomial(1, p[i])
    end
end

# ╔═╡ 7750a022-ef28-4104-84ac-ba43c03efaf5
begin
	m11_5_chain = sample(m11_5(chimpanzees.actor, chimpanzees.side, chimpanzees.cond, chimpanzees.pulled_left), NUTS(), 1000)
	m11_5_df = DataFrame(m11_5_chain)
	describe(m11_5_df)
end

# ╔═╡ ca353ff4-fabb-4c72-8ba4-4294def7a751
md" #### Code 11.20"

# ╔═╡ 3f4c5f34-216d-4b8e-84a1-2c6c8b83716d
l_fun1 = (r, (ai, bi, pull_left)) -> begin
    p = logistic(get(r, "a[$ai]", 0) + get(r, "b[$bi]", 0))
    binomlogpdf(1, p, pull_left)
end

# ╔═╡ f7940bb6-00f5-4ea3-87a2-cdaf719f6147
m11_4_ll_1 = link(m11_4_df, l_fun1, zip(chimpanzees.actor, chimpanzees.treatment, chimpanzees.pulled_left))

# ╔═╡ 9d3de9ec-3b74-41fc-86cd-4fca5ed0d3fe
m11_4_ll_2 = hcat(m11_4_ll_1...);

# ╔═╡ a734af10-a0ff-49bf-87ab-a39c12ba68a2
l_fun2 = (r, (ai, si, ci, pull_left)) -> begin
    p = logistic(get(r, "a[$ai]", 0) + get(r, "bs[$si]", 0) + get(r, "bc[$ci]", 0))
    binomlogpdf(1, p, pull_left)
end

# ╔═╡ 2b6cf6d6-5d97-46f5-8359-c069342afa82
# Code 11.21 and 11.22 were omitted, as they are stan-specific

# ╔═╡ 19a1560d-3524-484d-a4ca-e1f9678f3f4e
md" #### Code 11.23"

# ╔═╡ 1f67c25c-1ef2-4d07-aa04-e5e8e06b6092
mean(@. exp(m11_4_df."b[4]" - m11_4_df."b[2]"))

# ╔═╡ c70321cd-ef51-41bd-9744-aeec09bf2efb
md" #### Code 11.24"

# ╔═╡ 297964d6-d5ff-46a1-8b96-c6bba81f6be4
let
	gb = groupby(chimpanzees, [:treatment, :actor, :side, :cond])
	global d_aggregated = combine(gb, :pulled_left => sum => :left_pulls)
	first(d_aggregated, 8)
end

# ╔═╡ 5ce38815-7ef8-44fc-899e-e2158d6ca291
md" #### Code 11.25"

# ╔═╡ cbd9fc1e-6035-4513-b113-e36dd36c56b9
@model function m11_6(actor, treatment, left_pulls)
    act_levels = length(levels(actor))
    a ~ MvNormal(zeros(act_levels), 1.5)
    treat_levels = length(levels(treatment))
    b ~ MvNormal(zeros(treat_levels), 0.5)
    
    p = @. logistic(a[actor] + b[treatment])
    for i ∈ eachindex(left_pulls)
        left_pulls[i] ~ Binomial(18, p[i])
    end
end

# ╔═╡ a952f503-d7d7-4b31-8a5a-510770d15bea
begin
	m11_6_chain = sample(m11_6(d_aggregated.actor, d_aggregated.treatment, d_aggregated.left_pulls), NUTS(), 1000)
	m11_6_df = DataFrame(m11_6_chain)
end

# ╔═╡ 65c3888c-490f-46dd-9a41-c9fa135c5f3c
md" #### Code 11.26"

# ╔═╡ 216bb14d-f167-4b18-87e2-ae635a413622
l_fun = (r, (ai, bi, left_pulls)) -> begin
    p = logistic(get(r, "a[$ai]", 0) + get(r, "b[$bi]", 0))
    binomlogpdf(18, p, left_pulls)
end

# ╔═╡ 6e5e457f-e829-4bf4-95ee-cc71f5f05818
begin
	m11_5_ll_3 = link(m11_5_df, l_fun, zip(chimpanzees.actor, chimpanzees.side, chimpanzees.cond, chimpanzees.pulled_left))
	m11_5_ll_4 = hcat(m11_5_ll_3...);
	
	compare([m11_5_ll_4, m11_4_ll_2], :psis, mnames=["m5", "m4"])
end

# ╔═╡ a021c167-1334-4aaf-89b6-1b363ea90978
begin
	m11_6_ll = link(m11_6_df, l_fun, zip(d_aggregated.actor, d_aggregated.treatment, d_aggregated.left_pulls))
	m11_6_ll = hcat(m11_6_ll...)
	
	try
	    compare([m11_6_ll, m11_4_ll], :psis, mnames=["m6", "m4"])
	catch e
	    println(e)
	end
end

# ╔═╡ 087b4dea-68af-4022-a35b-b2e09e68ca59
md" #### Code 11.27"

# ╔═╡ ae7e3091-3917-4bc7-80d6-2770c72f374f
(
    -2*binomlogpdf(9, 0.2, 6),
    -2*sum(binomlogpdf.(1, 0.2, [1,1,1,1,1,1,0,0,0]))
)

# ╔═╡ e2403387-7305-4c58-98f8-c9e60061a23f
md" #### Code 11.28"

# ╔═╡ 51f74f4e-3b24-49ef-8221-0c74de9a8fd4
ucbadmit = CSV.read(sr_datadir("UCBadmit.csv"), DataFrame)

# ╔═╡ d2ef4d6c-6c77-4b88-8443-dd3dea5f4fa4
md" #### Code 11.29"

# ╔═╡ bf0f80c3-9049-4d0f-82cb-c2561fcc2aae
begin
	dat = ucbadmit[!, [:admit, :applications]]
	dat.gid = @. ifelse(ucbadmit.gender == "male", 1, 2)
	dat
end

# ╔═╡ 0290f649-11df-4d3e-81af-626d196316f3


# ╔═╡ 8b3dc9ea-fdac-4080-b284-4a8d3cf16f28
@model function m11_7(admit, applications, gid)
    a ~ MvNormal([0, 0], 1.5)
    p = @. logistic(a[gid])
    for i ∈ eachindex(applications)
        admit[i] ~ Binomial(applications[i], p[i])
    end
end

# ╔═╡ f60af3f5-1048-49fb-a3a7-3f0b2b2884f9
begin
	m11_7_chain = sample(m11_7(dat.admit, dat.applications, dat.gid), NUTS(), 1000)
	m11_7_df = DataFrame(m11_7_chain)
	describe(m11_7_df)
end

# ╔═╡ 67f5cc12-e357-4104-a3c6-7d3e381c281c
md" #### Code 11.30"

# ╔═╡ 540dcd1c-a44a-4c8a-9048-3cb5c382c278
let
	diff_a = m11_7_df."a[1]" .- m11_7_df."a[2]"
	diff_p = @. logistic(m11_7_df."a[1]") - logistic(m11_7_df."a[2]")
	describe(DataFrame(diff_a=diff_a, diff_p=diff_p))
end

# ╔═╡ 515ddb6a-5917-43c6-8312-96e57753b478
md" #### Code 11.31"

# ╔═╡ 151de157-f083-476b-8fb4-4f4814f0a062
let
	Random.seed!(1)
	
	fun = (r, (apps, gid)) -> begin
	    p = logistic(get(r, "a[$gid]", 0))
	    rand(Binomial(apps, p))
	end
	
	admit_pred = link(m11_7_df, fun, zip(dat.applications, dat.gid))
	admit_rate = @. admit_pred ./ dat.applications
	
	# plot
	xrange = 1:12
	p = plot(yrange=(0, 1), title="Posterior validation check", xlab="case", ylab="admit")
	scatter!(xrange, mean.(admit_rate), yerr=std.(admit_rate))
	pi_rate = PI.(admit_rate)
	scatter!(xrange, first.(pi_rate), shape=:cross, c=:black)
	scatter!(xrange, last.(pi_rate), shape=:cross, c=:black)
	
	dat_rates = dat.admit ./ dat.applications
	
	for idx in 1:2:12
	    r = [dat_rates[idx], dat_rates[idx+1]]
	    plot!([idx, idx+1], r, mark=:o, lw=2, c=:blue)
	    annotate!([(idx+0.5, mean(r)+0.03, (ucbadmit.dept[idx], 8))])
	end
	p
end

# ╔═╡ 22c6ef5a-6a54-4be6-a3ff-3a13382d5bcf
md" #### Code 11.32"

# ╔═╡ ef44deb9-b1f0-40a0-82a7-b970ed7fdfc0
dat.dept_id = reshape(hcat(1:6, 1:6)', 12)

# ╔═╡ 3ad96e80-177a-4b74-8d2f-9e1efc60ed4c
dat

# ╔═╡ 7bb690a7-04dc-4f0e-aabb-13190ec2cc6a
@model function model_m11_8(admit, applications, gid, dept_id)
    a ~ MvNormal([0, 0], 1.5)
    delta ~ MvNormal(zeros(6), 1.5)
    p = @. logistic(a[gid] + delta[dept_id])
    for i ∈ eachindex(applications)
        admit[i] ~ Binomial(applications[i], p[i])
    end
end

# ╔═╡ b9e52a64-ed27-49f0-b5d4-8f1d5af96457
begin
	m11_8n = sample(model_m11_8(dat.admit, dat.applications, dat.gid, dat.dept_id), NUTS(), 1000)
	m11_8_df = DataFrame(m11_8n)
	describe(m11_8_df)
end

# ╔═╡ ab439d9e-d427-4262-bacc-4e88759a04f4
md" #### Code 11.33"

# ╔═╡ 0e79c84d-c016-4ecd-bea1-b4a80f84e007
let
	diff_a = m11_8_df."a[1]" .- m11_8_df."a[2]"
	diff_p = logistic.(m11_8_df."a[1]") .- logistic.(m11_8_df."a[2]")
	describe(DataFrame(diff_a=diff_a, diff_p=diff_p))
end

# ╔═╡ cd0a1208-2721-4629-9784-1b1e137cd966
md" #### Code 11.34"

# ╔═╡ f768fbac-9982-44b8-9518-883760201de7
let
	gb1 = groupby(ucbadmit, :dept)
	fun = (as,gs) -> [(gender=g, a_ratio=a/sum(as)) for (a,g) in zip(as,gs)]
	c = combine(gb1, [:applications, :gender] => fun => AsTable)
	unstack(c, :gender, :dept, :a_ratio)
end

# ╔═╡ bb371f48-2aaf-4205-a3ae-832da7254c16
md" #### Code 11.32"

# ╔═╡ 989be969-0dd0-4ef9-9203-df59af4b7825
dat.dept_id = reshape(hcat(1:6, 1:6)', 12)

# ╔═╡ 46c1aa71-6177-4d9c-9fab-996cfe4187b9
@model function m11_8(admit, applications, gid, dept_id)
    a ~ MvNormal([0, 0], 1.5)
    delta ~ MvNormal(zeros(6), 1.5)
    p = @. logistic(a[gid] + delta[dept_id])
    for i ∈ eachindex(applications)
        admit[i] ~ Binomial(applications[i], p[i])
    end
end

# ╔═╡ 239721f2-2dcf-478c-9266-ff186389185f
let
	m11_8_chain = sample(m11_8(dat.admit, dat.applications, dat.gid, dat.dept_id), NUTS(), 1000)
	m11_8_df = DataFrame(m11_8_chain)
	describe(m11_8_df)
end

# ╔═╡ c0e44ae1-e1c9-481e-a748-e94d4fa4932b
md" #### Code 11.33"

# ╔═╡ 3bc11427-36c1-4a75-b8dd-36a61a7e8350
begin
	diff_a = m11_8_df."a[1]" .- m11_8_df."a[2]"
	diff_p = logistic.(m11_8_df."a[1]") .- logistic.(m11_8_df."a[2]")
	describe(DataFrame(diff_a=diff_a, diff_p=diff_p))
end

# ╔═╡ 9b4c1ba4-74c3-4ed0-a558-f326940335ff
# Code 11.34

gb = groupby(ucbadmit, :dept)

# ╔═╡ 7785e952-0275-44fd-a91a-fe08ee29368c
fun = (as,gs) -> [(gender=g, a_ratio=a/sum(as)) for (a,g) in zip(as,gs)]

# ╔═╡ baa3ad4e-539a-4f2a-9daa-c0183472ca66
c = combine(gb, [:applications, :gender] => fun => AsTable)

# ╔═╡ cd0358bb-f42e-4cd5-9dee-d0104cf3e07b
unstack(c, :gender, :dept, :a_ratio)

# ╔═╡ 0cafb1cb-956f-4ba2-afb6-3546c8f2d990
md" ## 11.2 Poisson regression."

# ╔═╡ 511ee325-6aab-4cb3-8e55-33f838007758
md" #### Code 11.35"

# ╔═╡ 766883c2-ff90-4ca3-bab4-278942bdf308
let
	Random.seed!(1)
	y_bin = rand(Binomial(1000, 1/1000), 10^5)
	mean(y_bin), var(y_bin)
end

# ╔═╡ 94195529-1d37-42cd-a5ee-9cf777b1f1bb
md" #### Code 11.36"

# ╔═╡ 8b4d82a3-8a9c-462f-a3cb-e5a71f94ec58
begin
	kline = CSV.read(sr_datadir("Kline.csv"), DataFrame)
	kline.P = standardize(ZScoreTransform, log.(kline.population))
	kline.contact_id = ifelse.(kline.contact .== "high", 2, 1)
end

# ╔═╡ b5d03baa-7e2c-403b-8d65-983dad4723f3
md" #### Code 11.38"

# ╔═╡ a49dc7d3-4743-4464-aabc-7a7fceb74bf6
let
	p = plot(xlims=(-1, 100), ylims=(0, 0.08), 
        xlab="mean number of tools", ylab="Density")
	r = 0:0.5:100
	plot!(r, LogNormal(0, 10), label=L"a \sim \ln \mathcal{N}(0, 10)")
	plot!(r, LogNormal(3, 0.5), label=L"a \sim \ln \mathcal{N}(3, 0.5)")
end

# ╔═╡ db5ad34c-65c4-432d-a6e2-c420d43d7faa
md" #### Code 11.39"

# ╔═╡ 8e16a5cb-1cd3-4ad5-8297-af28e503aac8
let
	a = rand(Normal(0, 10), 10^4)
	λ = exp.(a)
	mean(λ)
end

# ╔═╡ 6d083875-75e2-47b9-b58f-6cb6b5d95172
md" #### Code 11.40 - joined with 11.38"

# ╔═╡ 53acdf8b-a783-4272-bd5e-147a5b0e1ae8
md" #### Code 11.41"

# ╔═╡ 957254e4-c9c0-4954-b6c4-24e6146151ff
let
	Random.seed!(1)
	N = 100
	a = rand(Normal(3, 0.5), N)
	b = rand(Normal(0, 10), N)
	p = plot(xlims=(-2, 2), ylims=(0, 100))
	r = -2:0.05:2
	for (aᵢ, bᵢ) ∈ zip(a, b)
    	plot!(r, x -> exp(aᵢ + bᵢ * x), c=:gray)
	end
	p
end

# ╔═╡ 6c57ae79-8980-426a-bc27-d18b92e094e2
md" #### Code 11.42"

# ╔═╡ 6e4d7d10-1dae-40f6-9c1b-8ee0924c2f16
let
	Random.seed!(1)
	N = 100
	a = rand(Normal(3, 0.5), N)
	b = rand(Normal(0, .2), N)
	p = plot(xlims=(-2, 2), ylims=(0, 100))
	r = -2:0.05:2
	for (aᵢ, bᵢ) ∈ zip(a, b)
    	plot!(r, x -> exp(aᵢ + bᵢ * x), c=:gray)
	end
	p
end

# ╔═╡ 8d45b2c3-8305-477d-9f88-9c4152d726df
md" #### Code 11.43"

# ╔═╡ 99cee684-134c-4def-bb8c-c7cb03cd06ab
let
	N = 100
	x_seq = range(log(100), log(200_000), length=100)
	a = rand(Normal(3, 0.5), N)
	b = rand(Normal(0, .2), N)
	λ = map(x -> (@. exp(a + b * x)), x_seq)
	λ = hcat(λ...)
	p = plot(xlims=extrema(x_seq), ylims=(0, 500), xlab="log population", ylab="total tools")
	for λᵢ ∈ eachrow(λ)
    	plot!(x_seq, λᵢ, c=:gray)
	end
	p
end

# ╔═╡ 13c4376e-616f-44dc-8185-87a4f2def727
md" #### Code 11.44"

# ╔═╡ a5b54f8c-e670-4402-b218-8164170db469
let
	N = 100
	x_seq = range(log(100), log(200_000), length=100)
	a = rand(Normal(3, 0.5), N)
	b = rand(Normal(0, .2), N)
	λ = map(x -> (@. exp(a + b * x)), x_seq)
	x_seq = range(log(100), log(200_000), length=100)
	p = plot(xlims=extrema(exp.(x_seq)), ylims=(0, 500), 
	    xlab="population", ylab="total tools")
	for λᵢ ∈ eachrow(λ)
    	plot!(exp.(x_seq), λᵢ, c=:gray)
	end
	p
end

# ╔═╡ d0af8d51-7fc8-41cd-9a51-20775f11e1e1
md" #### Code 11.45"

# ╔═╡ 7d462c37-f172-4624-8012-02a44bb19a4e
# intercept only
@model function m11_9(T)
    a ~ Normal(3, 0.5)
    λ = exp(a)
    T ~ Poisson(λ)
end

# ╔═╡ 9e910804-5539-424f-a769-d22302020888
begin
	m11_9_ch = sample(m11_9(kline.total_tools), NUTS(), 1000)
	m11_9_df = DataFrame(m11_9_ch)
end

# ╔═╡ f603a8ff-ba98-4bc0-80c6-a89243f83a3d
# interaction model
@model function m11_10(T, cid, P)
    a ~ MvNormal([3, 3], 0.5)
    b ~ MvNormal([0, 0], 0.2)
    λ = @. exp(a[cid] + b[cid]*P)
    for i ∈ eachindex(T)
        T[i] ~ Poisson(λ[i])
    end
end

# ╔═╡ 36a6c6f8-3189-4cbe-a1de-dc93d02ccd51
m11_10_ch = sample(m11_10(kline.total_tools, kline.contact_id, kline.P), NUTS(), 1000)

# ╔═╡ 0646f466-3247-495a-8a7c-742256f64b11
m11_10_df = DataFrame(m11_10_ch);

# ╔═╡ 199010f9-bce8-4318-9912-3ed575ef21fe
md" #### Code 11.46"

# ╔═╡ 448b3506-728b-4631-b6cd-62df0cd39f45
let
	f = (r, T) -> poislogpdf(exp(r.a), T)
	ll_m9 = link(m11_9_df, f, kline.total_tools)
	ll_m9 = hcat(ll_m9...)
	f = (r, (T, cid, P)) -> poislogpdf(exp(get(r, "a[$cid]", 0) + get(r, "b[$cid]", 0)*P), T)
	ll_m10 = link(m11_10_df, f, zip(kline.total_tools, kline.contact_id, kline.P))
	global ll_m10 = hcat(ll_m10...)
	compare([ll_m9, ll_m10], :psis, mnames=["m9", "m10"])
end

# ╔═╡ 88c22cbd-5f7b-4e50-9a0e-0ab22b857559
md" #### Code 11.47"

# ╔═╡ 2b7adad6-be9f-45f8-9385-e5d31ed69067
let
	# get psis K values
	t = ll_m10'
	m10_t = collect(reshape(t, size(t)..., 1))
	PSIS_m10 = psis_loo(m10_t)
	k = PSIS_m10.pointwise(:pareto_k)
	mark_size = standardize(ZScoreTransform, k).+5
	mark_color = ifelse.(kline.contact_id .== 1, :white, :blue)
	scatter(kline.P, kline.total_tools, xlab="log population (std)", ylab="total tools", ylim=(0, 75), 
    	ms=mark_size, mc=mark_color, msw=2)
	ns = 100
	P_seq = range(-1.4, 3, length=ns)
	# cid=1 predictions
	f = (r,p) -> exp(r."a[1]" + r."b[1]"*p)
	λ = link(m11_10_df, f, P_seq)
	λ = hcat(λ...)
	λ1_mean = mean.(eachcol(λ))
	λ1_ci = PI.(eachcol(λ))
	λ1_ci = vcat(λ1_ci'...)
	plot!(P_seq, [λ1_mean λ1_mean], fillrange=λ1_ci, fillalpha=0.2, ls=:dash, c=:black, lw=1.5)
	f = (r,p) -> exp(r."a[2]" + r."b[2]"*p)
	# cid=2 predictions
	λ = link(m11_10_df, f, P_seq)
	λ = hcat(λ...)
	λ2_mean = mean.(eachcol(λ))
	λ2_ci = PI.(eachcol(λ))
	λ2_ci = vcat(λ2_ci'...)
	plot!(P_seq, [λ2_mean λ2_mean], fillrange=λ2_ci, fillalpha=0.2, c=:black, lw=1.5)
	# Code 11.48
	scatter(kline.population, kline.total_tools, xlab="population", ylab="total tools", ylim=(0, 75), 
    	ms=mark_size, mc=mark_color, msw=2)
	P_seq = range(-5, 3, length=ns)
	# 1.53 is sd of log(population)
	# 9 is mean of log(population)
	pop_seq = @. exp(P_seq * 1.53 + 9)
	plot!(pop_seq, [λ1_mean λ1_mean], fillrange=λ1_ci, fillalpha=0.2, ls=:dash, c=:black, lw=1.5)
	plot!(pop_seq, [λ2_mean λ2_mean], fillrange=λ2_ci, fillalpha=0.2, c=:black, lw=1.5)
end

# ╔═╡ fddece02-6e0e-466c-97aa-bb96fa1fe013
md" #### Code 11.49"

# ╔═╡ aff8e056-ba12-4e13-b6bf-63da9181b8b7
@model function m11_11(T, P, cid)
    a ~ MvNormal([1,1], 1)
    b₁ ~ Exponential(1)
    b₂ ~ Exponential(1)
    b = [b₁, b₂]
    g ~ Exponential(1)
    λ = @. exp(a[cid]) * P ^ b[cid] / g
    for i ∈ eachindex(T)
        T[i] ~ Poisson(λ[i])
    end
end

# ╔═╡ 421f0b89-bf68-4546-8184-445b2e2de9fe
m11_11_ch = sample(m11_11(kline.total_tools, kline.population, kline.contact_id), NUTS(), 1000)

# ╔═╡ 88648d82-1c06-4b8c-bc91-58518b707905
m11_11_df = DataFrame(m11_11_ch)

# ╔═╡ 3ab8e039-da2e-4cca-9a1c-94d929a7aa78
describe(m11_11_df)

# ╔═╡ ba34af6e-85ee-4ff1-95c1-f627405262a2
md" #### Code 11.50"

# ╔═╡ 4f6e6857-f4c7-4c28-94e5-d79390127cad
begin
	Random.seed!(1)
	num_days = 30
	y = rand(Poisson(1.5), num_days)
end

# ╔═╡ 6cb0c0a9-de92-44f8-8b76-931b8d7a8058
md" #### Code 11.51"

# ╔═╡ c87d89c9-f559-4210-a300-cd42ffdc0779
begin
	Random.seed!(2)
	num_weeks = 4
	y_new = rand(Poisson(0.5*7), num_weeks)
	y_all = [y; y_new]
end

# ╔═╡ 85909da0-d2b9-437e-ae47-44940fcd74df
md" #### Code 11.52"

# ╔═╡ ac56b077-a1c9-4e6d-9c78-64a613df5bb0
begin
	exposure = [repeat([1], 30); repeat([7], 4)]
	monastery = [repeat([0], 30); repeat([1], 4)]
	expmon = DataFrame(y=y_all, days=exposure, monastery=monastery)
end

# ╔═╡ 8d0e632f-03ea-439c-9f6e-d57f489fc92c
@model function m11_12(y, log_days, monastery)
    a ~ Normal()
    b ~ Normal()
    λ = @. exp(log_days + a + b*monastery)
    @. y ~ Poisson(λ)
end

# ╔═╡ 0c1ed9fa-80b0-449a-83be-f9e208f5e536
md" #### Code 11.53"

# ╔═╡ 545a7bd1-e56c-45e9-bdd6-efd3e14632fc
begin
	Random.seed!(1)
	expmon.log_days = log.(expmon.days)
	m11_12_chain = sample(m11_12(expmon.y, expmon.log_days, expmon.monastery), NUTS(), 1000)
	m11_12_df = DataFrame(m11_12_chain)
	describe(m11_12_df)
end

# ╔═╡ d83b7902-85fd-4424-8ffe-98f9da71d0ae
md" #### Code 11.54"

# ╔═╡ d4fa0864-c56d-4665-89a6-e9197e037c56
md"

!!! note
Values are slightly different from the book. This is due to non-Normal distributions (you can check this yourself with `optimize`)"

# ╔═╡ 65ed2cfe-4404-4756-a9d5-3d55cae1eefd
let
	λ_old = exp.(m11_12_df.a)
	λ_new = @. exp(m11_12_df.a + m11_12_df.b)
	describe(DataFrame(λ_old=λ_old, λ_new=λ_new))
end

# ╔═╡ 9d721d0d-abe9-4f42-bc98-6245aeb2b2d0
md" ## 11.3 Multinomial and categorical models"

# ╔═╡ 028d2911-1c4a-4de1-ae4d-83e8fb153f12
md" #### Code 11.55"

# ╔═╡ a9c9cc35-05f6-49c8-a3df-2f303cf7885c
md" #### Code 11.56"

# ╔═╡ bc2d4336-7a19-4179-b791-e0cdb63cee0a
@model function m11_13(career, income)
    a ~ MvNormal([0, 0], 1)
    b ~ TruncatedNormal(0, 0.5, 0, Inf)
    p = softmax([a[1] + b*income[1], a[2] + b*income[2], 0])
    career ~ Categorical(p)
end

# ╔═╡ ac7531ca-02b5-44e9-a418-f8c399e668e8
md" #### Code 11.57"

# ╔═╡ 7690e5a3-0473-48fa-963c-7b55d00a4105
let
	# simulate career choice among 500 individuals
	N = 500
	global income = [1,2,5]
	c_score = 0.5*income
	p = softmax(c_score)
	w = Weights(p)
	
	# simulate choice
	Random.seed!(34302)
	career = [sample(w) for _ in 1:N]
	Random.seed!(121)
	m11_13_chain = sample(m11_13(career, income), NUTS(), 5000, n_chains=4)
	global m11_13_df = DataFrame(m11_13_chain)
	describe(m11_13_df)
end

# ╔═╡ 37e13f5d-0d82-4c2d-a1c2-479843eb7886
md" #### Code 11.58"

# ╔═╡ 3e8dd9ee-867e-43e5-b74c-d4f586bfd41b
let
	# logit scores
	s1 = m11_13_df."a[1]" + m11_13_df.b * income[1]
	s2_orig = m11_13_df."a[2]" + m11_13_df.b * income[2]
	s2_new = m11_13_df."a[2]" + m11_13_df.b * income[2] * 2
	# compute probabilities for original and counterfactual
	p_orig = softmax.(eachrow(hcat(s1, s2_orig, zeros(length(s1)))))
	p_orig = hcat(p_orig...)'
	p_new = softmax.(eachrow(hcat(s1, s2_new, zeros(length(s1)))))
	p_new = hcat(p_new...)'
	p_diff = p_new[:, 2] - p_orig[:, 2]
	describe(DataFrame(p_diff=p_diff))
end

# ╔═╡ ed01aaf1-1c1e-4441-9999-1a04916df98e
@model function m11_14(career, family_income)
    a ~ MvNormal([0, 0], 1.5)
    b ~ MvNormal([0, 0], 1)
    
    for i ∈ eachindex(career)
        income = family_income[i]
        s = [a[1] + b[1] * income, a[2] + b[2] * income, 0]
        p = softmax(s)
        career[i] ~ Categorical(p)
    end
end

# ╔═╡ 13d724a6-6902-4ec1-86e8-3c6a86e3957d
md" #### Code 11.59"

# ╔═╡ aff70993-542e-4be7-99cb-160059a84e06
let
	Random.seed!(1)
	N = 500
	family_income = rand(Uniform(), N)
	b = [-2, 0, 2]
	career = [
    begin
        score = (0.5 * 1:3) + b * inc
        p = softmax(score)
        sample(Weights(p))
    end
    for inc in family_income
	]
	m11_14_chain = sample(m11_14(career, family_income), NUTS(), 1000)
	m11_14_df = DataFrame(m11_14_chain)
	describe(m11_14_df)
end

# ╔═╡ dae83135-14da-4b46-b974-782d4ec4e4d9
md" #### Code 11.61"

# ╔═╡ 4f2e5271-e301-4f44-8474-ec9b5c97150c
@model function m_binom(admit, applications)
    a ~ Normal(0, 1.5)
    p = logistic(a)
    @. admit ~ Binomial(applications, p)
end

# ╔═╡ 2f08dbb9-f31a-4fc2-ab11-fdfb88e5ddb9
m_binom_chain = sample(m_binom(ucbadmit.admit, ucbadmit.applications), NUTS(), 1000)

# ╔═╡ cc4700af-9c54-4060-a9ce-8c426e25eeb1
m_binom_df = DataFrame(m_binom_chain);

# ╔═╡ cf9ade41-251c-4320-a606-4e877b996632
@model function m_pois(admit, reject)
    a ~ MvNormal([0, 0], 1.5)
    λ = exp.(a)
    admit ~ Poisson(λ[1])
    reject ~ Poisson(λ[2])
end

# ╔═╡ d6f7eff0-9142-4635-81f3-11ec8c2a500b
m_pois_chain = sample(m_pois(ucbadmit.admit, ucbadmit.reject), NUTS(), 1000)

# ╔═╡ 06cf33fc-f0f8-4cc7-b811-d035fbf40b20
m_pois_df = DataFrame(m_pois_chain);

# ╔═╡ d6fa7523-c369-454d-bce9-a12132fe01f1
md" #### Code 11.62"

# ╔═╡ e2fbc671-9ea3-4ba3-b847-b71ed7a67933
m_binom_df.a .|> logistic |> mean

# ╔═╡ 17745763-875a-45f3-9b1a-8c4c80d5060c
md" #### Code 11.63"

# ╔═╡ eca35612-759b-43ef-a94b-69dc4856317f
a1 = m_pois_df."a[1]" |> mean

# ╔═╡ 040be4fa-9158-45da-a9cd-ac419d59bb56
a2 = m_pois_df."a[2]" |> mean

# ╔═╡ 4b34baa2-9528-4311-a54d-33f080f64103
exp(a1)/(exp(a1)+exp(a2))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DrWatson = "634d3b9d-ee7a-5ddf-bec9-22491ea816e1"
FreqTables = "da1fdf0e-e0ff-5433-a45f-9bb5ff651cb1"
Logging = "56ddb016-857b-54e1-b83d-db4d58db5568"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
ParetoSmooth = "a68b5a21-f429-434e-8bfa-46b447300aac"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatisticalRethinking = "2d09df54-9d0f-5258-8220-54c2a3d4fbee"
StatisticalRethinkingPlots = "e1a513d0-d9d9-49ff-a6dd-9d2e9db473da"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"
Turing = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"

[compat]
CSV = "~0.10.12"
DataFrames = "~1.6.1"
DrWatson = "~2.13.0"
FreqTables = "~0.4.6"
Optim = "~1.8.0"
ParetoSmooth = "~0.7.4"
StatisticalRethinking = "~4.7.0"
StatisticalRethinkingPlots = "~1.1.0"
StatsBase = "~0.33.21"
StatsPlots = "~0.15.6"
Turing = "~0.28.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "5ba12210c7c993090f53a79d9741b99495374bb0"

[[deps.ADTypes]]
git-tree-sha1 = "41c37aa88889c171f1300ceac1313c06e891d245"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.6"

[[deps.ANSIColoredPrinters]]
git-tree-sha1 = "574baf8110975760d391c710b6341da1afa48d8c"
uuid = "a4c015fc-c6ff-483c-b24f-f7ea428134e9"
version = "0.0.1"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractMCMC]]
deps = ["BangBang", "ConsoleProgressMonitor", "Distributed", "LogDensityProblems", "Logging", "LoggingExtras", "ProgressLogging", "Random", "StatsBase", "TerminalLoggers", "Transducers"]
git-tree-sha1 = "87e63dcb990029346b091b170252f3c416568afc"
uuid = "80f14c24-f653-4e6a-9b94-39d6b0f70001"
version = "4.4.2"

[[deps.AbstractPPL]]
deps = ["AbstractMCMC", "DensityInterface", "Random", "Setfield", "SparseArrays"]
git-tree-sha1 = "caa9b62583577b0d6b222f11f54aa29fabbdb5ca"
uuid = "7a57a42e-76ec-4ea3-a279-07e840d6d9cf"
version = "0.6.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Test"]
git-tree-sha1 = "cb96992f1bec110ad211b7e410e57ddf7944c16f"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.35"
weakdeps = ["AxisKeys", "IntervalSets", "Requires", "StaticArrays", "StructArrays"]

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cde29ddf7e5726c9fb511f340244ea3481267608"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdvancedHMC]]
deps = ["AbstractMCMC", "ArgCheck", "DocStringExtensions", "InplaceOps", "LinearAlgebra", "LogDensityProblems", "LogDensityProblemsAD", "ProgressMeter", "Random", "Requires", "Setfield", "SimpleUnPack", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "acbe805c3078ba0057bb56985248bd66bce016b1"
uuid = "0bf59076-c3b1-5ca4-86bd-e02cd72cde3d"
version = "0.5.5"

    [deps.AdvancedHMC.extensions]
    AdvancedHMCCUDAExt = "CUDA"
    AdvancedHMCMCMCChainsExt = "MCMCChains"
    AdvancedHMCOrdinaryDiffEqExt = "OrdinaryDiffEq"

    [deps.AdvancedHMC.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
    OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"

[[deps.AdvancedMH]]
deps = ["AbstractMCMC", "Distributions", "FillArrays", "LinearAlgebra", "LogDensityProblems", "Random", "Requires"]
git-tree-sha1 = "b2a1602952739e589cf5e2daff1274a49f22c9a4"
uuid = "5b7e9947-ddc0-4b3f-9b55-0d8042f74170"
version = "0.7.5"
weakdeps = ["DiffResults", "ForwardDiff", "MCMCChains", "StructArrays"]

    [deps.AdvancedMH.extensions]
    AdvancedMHForwardDiffExt = ["DiffResults", "ForwardDiff"]
    AdvancedMHMCMCChainsExt = "MCMCChains"
    AdvancedMHStructArraysExt = "StructArrays"

[[deps.AdvancedPS]]
deps = ["AbstractMCMC", "Distributions", "Libtask", "Random", "Random123", "StatsFuns"]
git-tree-sha1 = "4d73400b3583147b1b639794696c78202a226584"
uuid = "576499cb-2369-40b2-a588-c64705576edc"
version = "0.4.3"

[[deps.AdvancedVI]]
deps = ["Bijectors", "Distributions", "DistributionsAD", "DocStringExtensions", "ForwardDiff", "LinearAlgebra", "ProgressMeter", "Random", "Requires", "StatsBase", "StatsFuns", "Tracker"]
git-tree-sha1 = "1f919a9c59cf3dfc68b64c22c453a2e356fca473"
uuid = "b5ca4192-6429-45e5-a2d9-87aec30a685c"
version = "0.2.4"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "bbec08a37f8722786d87bedf84eae19c020c4efa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.7.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "c06a868224ecba914baa6942988e2f2aade419be"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "0.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.AxisKeys]]
deps = ["AbstractFFTs", "ChainRulesCore", "CovarianceEstimation", "IntervalSets", "InvertedIndices", "LazyStack", "LinearAlgebra", "NamedDims", "OffsetArrays", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "dba0fdaa3a95e591aa9cbe0df9aba41e295a2011"
uuid = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
version = "0.2.13"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables"]
git-tree-sha1 = "7aa7ad1682f3d5754e3491bb59b8103cae28e3a3"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.40"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "f1f03a9fa24271160ed7e73051fba3c1a759b53f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.4.0"

[[deps.Bijectors]]
deps = ["ArgCheck", "ChainRules", "ChainRulesCore", "ChangesOfVariables", "Compat", "Distributions", "Functors", "InverseFunctions", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "MappedArrays", "Random", "Reexport", "Requires", "Roots", "SparseArrays", "Statistics"]
git-tree-sha1 = "199dc2c4151db557549a0ad8888ce1a60337ff42"
uuid = "76274a88-744f-5084-9051-94815aaf08c4"
version = "0.13.8"

    [deps.Bijectors.extensions]
    BijectorsDistributionsADExt = "DistributionsAD"
    BijectorsForwardDiffExt = "ForwardDiff"
    BijectorsLazyArraysExt = "LazyArrays"
    BijectorsReverseDiffExt = "ReverseDiff"
    BijectorsTrackerExt = "Tracker"
    BijectorsZygoteExt = "Zygote"

    [deps.Bijectors.weakdeps]
    DistributionsAD = "ced4e74d-a319-5a8a-b0ac-84af2272839c"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "601f7e7b3d36f18790e2caf83a882d88e9b71ff1"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.4"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "679e69c611fff422038e9e21e270c4197d49d918"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.12"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "1568b28f91293458345dabba6a5ea3f183250a61"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.8"

    [deps.CategoricalArrays.extensions]
    CategoricalArraysJSONExt = "JSON"
    CategoricalArraysRecipesBaseExt = "RecipesBase"
    CategoricalArraysSentinelArraysExt = "SentinelArrays"
    CategoricalArraysStructTypesExt = "StructTypes"

    [deps.CategoricalArrays.weakdeps]
    JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SentinelArrays = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
    StructTypes = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"

[[deps.ChainRules]]
deps = ["Adapt", "ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "SparseInverseSubset", "Statistics", "StructArrays", "SuiteSparse"]
git-tree-sha1 = "4cfc4916725289132746f730755262e1919cff38"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.59.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "0d12ee16b3f62e4e33c3277773730a5b21a74152"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.20.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "2fba81a302a7be671aefe194f0525ef231104e7f"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.8"
weakdeps = ["InverseFunctions"]

    [deps.ChangesOfVariables.extensions]
    ChangesOfVariablesInverseFunctionsExt = "InverseFunctions"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "9ebb045901e9bbf58767a9f34ff89831ed711aae"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.7"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "c0ae2a86b162fb5d7acc65269b469ff5b8a73594"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.8.1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "cd67fc487743b2f0fd4380d4cbd3a24660d0eec8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.3"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "75bd5b6fc5089df449b5d35fa501c846c9b6549b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.12.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "8cfa272e8bdedfa88b6aefbbca7c19f1befac519"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.3.0"

[[deps.ConsoleProgressMonitor]]
deps = ["Logging", "ProgressMeter"]
git-tree-sha1 = "3ab7b2136722890b9af903859afcf457fa3059e8"
uuid = "88cd18e8-d9cc-4ea6-8889-5259c0d15c8b"
version = "0.1.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.CovarianceEstimation]]
deps = ["LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "9a44ddc9e60ee398934b73a5168f5806989e6792"
uuid = "587fd27a-f159-11e8-2dae-1979310e6154"
version = "0.2.11"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "ac67408d9ddf207de5cfa9a97e114352430f01ed"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.16"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "7c302d7a5fec5214eb8a5a4c466dcf7a51fcf169"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.107"
weakdeps = ["ChainRulesCore", "DensityInterface", "Test"]

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

[[deps.DistributionsAD]]
deps = ["Adapt", "ChainRules", "ChainRulesCore", "Compat", "Distributions", "FillArrays", "LinearAlgebra", "PDMats", "Random", "Requires", "SpecialFunctions", "StaticArrays", "StatsFuns", "ZygoteRules"]
git-tree-sha1 = "d61f08c7bd15c5ab215fd7a2eb61c1ae15d8ff5e"
uuid = "ced4e74d-a319-5a8a-b0ac-84af2272839c"
version = "0.6.53"

    [deps.DistributionsAD.extensions]
    DistributionsADForwardDiffExt = "ForwardDiff"
    DistributionsADLazyArraysExt = "LazyArrays"
    DistributionsADReverseDiffExt = "ReverseDiff"
    DistributionsADTrackerExt = "Tracker"

    [deps.DistributionsAD.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Documenter]]
deps = ["ANSIColoredPrinters", "Base64", "Dates", "DocStringExtensions", "IOCapture", "InteractiveUtils", "JSON", "LibGit2", "Logging", "Markdown", "REPL", "Test", "Unicode"]
git-tree-sha1 = "39fd748a73dce4c05a9655475e437170d8fb1b67"
uuid = "e30172f5-a6a5-5a46-863b-614d45cd2de4"
version = "0.27.25"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DrWatson]]
deps = ["Dates", "FileIO", "JLD2", "LibGit2", "MacroTools", "Pkg", "Random", "Requires", "Scratch", "UnPack"]
git-tree-sha1 = "f83dbe0ef99f1cf32b815f0dad632cb25129604e"
uuid = "634d3b9d-ee7a-5ddf-bec9-22491ea816e1"
version = "2.13.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPPL]]
deps = ["AbstractMCMC", "AbstractPPL", "BangBang", "Bijectors", "ChainRulesCore", "ConstructionBase", "Distributions", "DocStringExtensions", "LinearAlgebra", "LogDensityProblems", "MacroTools", "OrderedCollections", "Random", "Setfield", "Test", "ZygoteRules"]
git-tree-sha1 = "a32e74c97f2e6e087a6870f7a641e84621d7e9e9"
uuid = "366bfd00-2699-11ea-058f-f148b4cae6d8"
version = "0.23.11"

[[deps.EllipticalSliceSampling]]
deps = ["AbstractMCMC", "ArrayInterface", "Distributions", "Random", "Statistics"]
git-tree-sha1 = "973b4927d112559dc737f55d6bf06503a5b3fc14"
uuid = "cad2338a-1db2-11e9-3401-43bc07c9ede2"
version = "1.1.0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "c5c28c245101bd59154f649e19b038d15901b5dc"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.2"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "9f00e42f8d99fdde64d40c8ea5d14269a2e2c1aa"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.21"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "5b93957f6dcd33fc343044af3d48c215be2562f1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.9.3"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "73d1214fec245096717847c62d389a5d2ac86504"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.22.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FreqTables]]
deps = ["CategoricalArrays", "Missings", "NamedArrays", "Tables"]
git-tree-sha1 = "4693424929b4ec7ad703d68912a6ad6eff103cfe"
uuid = "da1fdf0e-e0ff-5433-a45f-9bb5ff651cb1"
version = "0.4.6"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9a68d75d466ccc1218d0552a8e1631151c569545"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.5"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "a8c834cdae6a8347c72eea19930ebdaabb6015e6"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.1"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2abcce0c099dfb0863efc261be904fc2b85eccdd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.1+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "abbbb9ec3afd783a7cbd82ef01dcd088ea051398"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.1"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "eb8fed28f4994600e29beef49744639d985a04b2"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.16"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InplaceOps]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "50b41d59e7164ab6fda65e71049fee9d890731ff"
uuid = "505f98c9-085e-5b2c-8e89-488be7bf1f34"
version = "0.3.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5fdf2fe6724d8caabf43b557b84ce53f3b7e2f6b"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.0.2+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalSets]]
deps = ["Dates", "Random"]
git-tree-sha1 = "3d8866c029dd6b16e69e0d4a939c4dfcb98fac47"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.8"
weakdeps = ["Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "68772f49f54b479fa88ace904f6127f0a3bb2e46"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.12"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "PrecompileTools", "Printf", "Reexport", "Requires", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "7c0008f0b7622c6c0ee5c65cbc667b69f8a65672"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.45"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60b1194df0a3298f460063de985eae7b01bc011a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.1+0"

[[deps.KernelAbstractions]]
deps = ["Adapt", "Atomix", "InteractiveUtils", "LinearAlgebra", "MacroTools", "PrecompileTools", "Requires", "SparseArrays", "StaticArrays", "UUIDs", "UnsafeAtomics", "UnsafeAtomicsLLVM"]
git-tree-sha1 = "4e0cb2f5aad44dcfdc91088e85dee4ecb22c791c"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.16"

    [deps.KernelAbstractions.extensions]
    EnzymeExt = "EnzymeCore"

    [deps.KernelAbstractions.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "fee018a29b60733876eb557804b5b109dd3dd8a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.8"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Preferences", "Printf", "Requires", "Unicode"]
git-tree-sha1 = "cb4619f7353fc62a1a22ffa3d7ed9791cfb47ad8"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "6.4.2"

    [deps.LLVM.extensions]
    BFloat16sExt = "BFloat16s"

    [deps.LLVM.weakdeps]
    BFloat16s = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "98eaee04d96d973e79c25d49167668c5c8fb50e2"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.27+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LRUCache]]
git-tree-sha1 = "b3cc6698599b10e652832c2f23db3cab99d51b59"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.6.1"
weakdeps = ["Serialization"]

    [deps.LRUCache.extensions]
    SerializationExt = ["Serialization"]

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "62edfee3211981241b57ff1cedf4d74d79519277"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.15"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyStack]]
deps = ["ChainRulesCore", "LinearAlgebra", "NamedDims", "OffsetArrays"]
git-tree-sha1 = "2eb4a5bf2eb0519ebf40c797ba5637d327863637"
uuid = "1fad7336-0346-5a1a-a56f-a06ba010965b"
version = "0.0.8"

[[deps.LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "fb6803dafae4a5d62ea5cab204b1e657d9737e7f"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.2.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtask]]
deps = ["FunctionWrappers", "LRUCache", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "345a40c746404dd9cb1bbc368715856838ab96f2"
uuid = "6f1fad26-d15e-5dc8-ae53-837a1d7b8c9f"
version = "0.8.6"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogDensityProblems]]
deps = ["ArgCheck", "DocStringExtensions", "Random"]
git-tree-sha1 = "f9a11237204bc137617194d79d813069838fcf61"
uuid = "6fdf6af0-433a-55f7-b3ed-c6c6e0b8df7c"
version = "2.1.1"

[[deps.LogDensityProblemsAD]]
deps = ["DocStringExtensions", "LogDensityProblems", "Requires", "SimpleUnPack"]
git-tree-sha1 = "9c50732cd0f188766b6217ed6a2ebbdaf9890029"
uuid = "996a588d-648d-4e1f-a8f0-a84b347e47b1"
version = "1.7.0"

    [deps.LogDensityProblemsAD.extensions]
    LogDensityProblemsADADTypesExt = "ADTypes"
    LogDensityProblemsADEnzymeExt = "Enzyme"
    LogDensityProblemsADFiniteDifferencesExt = "FiniteDifferences"
    LogDensityProblemsADForwardDiffBenchmarkToolsExt = ["BenchmarkTools", "ForwardDiff"]
    LogDensityProblemsADForwardDiffExt = "ForwardDiff"
    LogDensityProblemsADReverseDiffExt = "ReverseDiff"
    LogDensityProblemsADTrackerExt = "Tracker"
    LogDensityProblemsADZygoteExt = "Zygote"

    [deps.LogDensityProblemsAD.weakdeps]
    ADTypes = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
    BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"
weakdeps = ["ChainRulesCore", "ChangesOfVariables", "InverseFunctions"]

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MCMCChains]]
deps = ["AbstractMCMC", "AxisArrays", "Compat", "Dates", "Distributions", "Formatting", "IteratorInterfaceExtensions", "KernelDensity", "LinearAlgebra", "MCMCDiagnosticTools", "MLJModelInterface", "NaturalSort", "OrderedCollections", "PrettyTables", "Random", "RecipesBase", "Serialization", "Statistics", "StatsBase", "StatsFuns", "TableTraits", "Tables"]
git-tree-sha1 = "f5f347b828fd95ece7398f412c81569789361697"
uuid = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
version = "5.5.0"

[[deps.MCMCDiagnosticTools]]
deps = ["AbstractFFTs", "DataAPI", "Distributions", "LinearAlgebra", "MLJModelInterface", "Random", "SpecialFunctions", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "dab8e9e9cd714fd3adc2344a97504e7e64e78546"
uuid = "be115224-59cd-429b-ad48-344e309966f0"
version = "0.1.5"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "72dc3cf284559eb8f53aa593fe62cb33f83ed0c0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.0.0+0"

[[deps.MLJModelInterface]]
deps = ["Random", "ScientificTypesBase", "StatisticalTraits"]
git-tree-sha1 = "14bd8088cf7cd1676aa83a57004f8d23d43cd81e"
uuid = "e80e1ace-859a-464e-9ed9-23947d8ae3ea"
version = "1.9.5"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "PrecompileTools", "Printf", "SparseArrays", "SpecialFunctions", "Test", "Unicode"]
git-tree-sha1 = "e2ae8cf5ac6daf5a3959f7f6ded9c2028b61d09d"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.25.1"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "629afd7d10dbc6935ec59b32daeb33bc4460a42e"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.4"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MonteCarloMeasurements]]
deps = ["Distributed", "Distributions", "ForwardDiff", "GenericSchur", "LinearAlgebra", "MacroTools", "Random", "RecipesBase", "Requires", "SLEEFPirates", "StaticArrays", "Statistics", "StatsBase", "Test"]
git-tree-sha1 = "19d4a73e20ca54f0f0e8a4ed349ee0dfd6e997b7"
uuid = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
version = "1.1.6"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "68bf5103e002c44adfd71fea6bd770b3f0586843"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.2"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "806eea990fb41f9b36f1253e5697aa645bf6a9f8"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.4.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NNlib]]
deps = ["Adapt", "Atomix", "ChainRulesCore", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "Pkg", "Random", "Requires", "Statistics"]
git-tree-sha1 = "900a11b3a2b02e36b25cb55a80777d4a4670f0f6"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.9.10"

    [deps.NNlib.extensions]
    NNlibAMDGPUExt = "AMDGPU"
    NNlibCUDACUDNNExt = ["CUDA", "cuDNN"]
    NNlibCUDAExt = "CUDA"
    NNlibEnzymeCoreExt = "EnzymeCore"

    [deps.NNlib.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    cuDNN = "02a925ec-e4fe-4b08-9a7e-0d78e3d38ccd"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NamedArrays]]
deps = ["Combinatorics", "DataStructures", "DelimitedFiles", "InvertedIndices", "LinearAlgebra", "Random", "Requires", "SparseArrays", "Statistics"]
git-tree-sha1 = "b84e17976a40cb2bfe3ae7edb3673a8c630d4f95"
uuid = "86f7a689-2022-50b4-a561-43c23ac3c673"
version = "0.9.8"

[[deps.NamedDims]]
deps = ["AbstractFFTs", "ChainRulesCore", "CovarianceEstimation", "LinearAlgebra", "Pkg", "Requires", "Statistics"]
git-tree-sha1 = "cb8ebcee2b4e07b72befb9def593baef8aa12f07"
uuid = "356022a1-0364-5f58-8944-0da4b18d706f"
version = "0.2.50"

[[deps.NamedTupleTools]]
git-tree-sha1 = "90914795fc59df44120fe3fff6742bb0d7adb1d0"
uuid = "d9ec5142-1e00-5aa0-9d6a-321866360f50"
version = "0.14.3"

[[deps.NaturalSort]]
git-tree-sha1 = "eda490d06b9f7c00752ee81cfa451efe55521e21"
uuid = "c020b1a1-e9b0-503a-9c33-f039bfc54a85"
version = "1.0.0"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ded64ff6d4fdd1cb68dfcbb818c69e144a5b2e4c"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.16"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "6a731f2b5c03157418a20c12195eb4b74c8f8621"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.13.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc6e1927ac521b659af340e0ca45828a3ffc748f"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.12+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "MathOptInterface", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "f55af9918e2a67dcadf5ec758a5ff25746c3819f"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.8.0"

[[deps.Optimisers]]
deps = ["ChainRulesCore", "Functors", "LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "34205b1204cc83c43cd9cfe53ffbd3b310f6e8c5"
uuid = "3bd65402-5787-11e9-1adc-39752487f4e2"
version = "0.3.1"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.ParetoSmooth]]
deps = ["AxisKeys", "InteractiveUtils", "LinearAlgebra", "LogExpFunctions", "MCMCDiagnosticTools", "NamedDims", "PrettyTables", "Printf", "Random", "Requires", "Statistics", "StatsBase"]
git-tree-sha1 = "fa3adfaaf140069ab48a2f1a3216fc5d45c06721"
uuid = "a68b5a21-f429-434e-8bfa-46b447300aac"
version = "0.7.4"

[[deps.ParetoSmoothedImportanceSampling]]
deps = ["CSV", "DataFrames", "Distributions", "JSON", "Printf", "Random", "Statistics", "StatsFuns", "Test"]
git-tree-sha1 = "c678e21715f9b6bbf4cc63047f935a68a9b44f20"
uuid = "98f080ec-61e2-11eb-1c7b-31ea1097256f"
version = "1.5.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "862942baf5663da528f66d24996eb6da85218e76"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "38a748946dca52a622e79eea6ed35c6737499109"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "88b895d13d53b5577fd53379d913b9ab9ac82660"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "00099623ffee15972c16111bcf84c58a0051257c"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.9.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "c860e84651f58ce240dd79e5d9e055d55234c35a"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.6.2"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "d7087c013e8a496ff396bae843b1e16d9a30ede8"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.10"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Roots]]
deps = ["Accessors", "ChainRulesCore", "CommonSolve", "Printf"]
git-tree-sha1 = "39ebae5b76c8cd5629bec21adfca78b437dac1e6"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.1.1"

    [deps.Roots.extensions]
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "6aacc5eefe8415f47b3e34214c1d79d2674a0ba2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.12"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "3aac6d68c5e57449f5b9b865c9ba50ac2970c4cf"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.42"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "ChainRulesCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FillArrays", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "916b8a94c0d61fa5f7c5295649d3746afb866aff"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.98.1"

    [deps.SciMLBase.extensions]
    ZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "51ae235ff058a64815e0a2c34b1db7578a06813d"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.7"

[[deps.ScientificTypesBase]]
git-tree-sha1 = "a8e18eb383b5ecf1b5e6fc237eb39255044fd92b"
uuid = "30f210dd-8aff-4c5f-94ba-8e64358c1161"
version = "3.0.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "0e7508ff27ba32f26cd459474ca2ede1bc10991f"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseInverseSubset]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "52962839426b75b3021296f7df242e40ecfc0852"
uuid = "dc90abb0-5640-4711-901d-7e5b23a2fada"
version = "0.1.2"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "f295e0a1da4ca425659c57441bcb59abb035a4bc"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.8"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Requires", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "5d66818a39bb04bf328e92bc933ec5b4ee88e436"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.5.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "f68dd04d131d9a8a8eb836173ee8f105c360b0c5"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.1"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.StatisticalRethinking]]
deps = ["CSV", "DataFrames", "Dates", "Distributions", "DocStringExtensions", "Documenter", "Formatting", "KernelDensity", "LinearAlgebra", "MCMCChains", "MonteCarloMeasurements", "NamedArrays", "NamedTupleTools", "Optim", "OrderedCollections", "Parameters", "ParetoSmoothedImportanceSampling", "PrettyTables", "Random", "Reexport", "Requires", "Statistics", "StatsBase", "StatsFuns", "StructuralCausalModels", "Tables", "Test", "Unicode"]
git-tree-sha1 = "0e4f7aff217dff83b0d2ddb883989dd01d23735a"
uuid = "2d09df54-9d0f-5258-8220-54c2a3d4fbee"
version = "4.7.0"

[[deps.StatisticalRethinkingPlots]]
deps = ["Distributions", "DocStringExtensions", "KernelDensity", "LaTeXStrings", "Parameters", "Plots", "Reexport", "Requires", "StatisticalRethinking", "StatsPlots"]
git-tree-sha1 = "bd7bd318815654491e6350c662020119be7792a0"
uuid = "e1a513d0-d9d9-49ff-a6dd-9d2e9db473da"
version = "1.1.0"

[[deps.StatisticalTraits]]
deps = ["ScientificTypesBase"]
git-tree-sha1 = "30b9236691858e13f167ce829490a68e1a597782"
uuid = "64bff920-2084-43da-a3e6-9bb72801c0c9"
version = "3.2.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "9115a29e6c2cf66cf213ccc17ffd61e27e743b24"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.6"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.StructArrays]]
deps = ["Adapt", "ConstructionBase", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "1b0b1205a56dc288b71b1961d48e351520702e24"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.17"

[[deps.StructuralCausalModels]]
deps = ["CSV", "Combinatorics", "DataFrames", "DataStructures", "Distributions", "DocStringExtensions", "LinearAlgebra", "NamedArrays", "Reexport", "Statistics"]
git-tree-sha1 = "01c838be8d7119708b839aa16d413088a1076ee8"
uuid = "a41e6734-49ce-4065-8b83-aff084c01dfd"
version = "1.4.1"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TerminalLoggers]]
deps = ["LeftChildRightSiblingTrees", "Logging", "Markdown", "Printf", "ProgressLogging", "UUIDs"]
git-tree-sha1 = "f133fab380933d042f6796eda4e130272ba520ca"
uuid = "5d786b92-1e48-4d6f-9151-6b4477ca9bed"
version = "0.1.7"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tracker]]
deps = ["Adapt", "DiffRules", "ForwardDiff", "Functors", "LinearAlgebra", "LogExpFunctions", "MacroTools", "NNlib", "NaNMath", "Optimisers", "Printf", "Random", "Requires", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "5c942be30a85ac75d14e9e527b55504031e1bbd3"
uuid = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
version = "0.2.31"
weakdeps = ["PDMats"]

    [deps.Tracker.extensions]
    TrackerPDMatsExt = "PDMats"

[[deps.TranscodingStreams]]
git-tree-sha1 = "1fbeaaca45801b4ba17c251dd8603ef24801dd84"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.2"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "3064e780dbb8a9296ebb3af8f440f787bb5332af"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.80"

    [deps.Transducers.extensions]
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.Turing]]
deps = ["AbstractMCMC", "AdvancedHMC", "AdvancedMH", "AdvancedPS", "AdvancedVI", "BangBang", "Bijectors", "DataStructures", "Distributions", "DistributionsAD", "DocStringExtensions", "DynamicPPL", "EllipticalSliceSampling", "ForwardDiff", "Libtask", "LinearAlgebra", "LogDensityProblems", "LogDensityProblemsAD", "MCMCChains", "NamedArrays", "Printf", "Random", "Reexport", "Requires", "SciMLBase", "Setfield", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Tracker"]
git-tree-sha1 = "5e9b655060490d1a16cbe569acb276578ace48aa"
uuid = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"
version = "0.28.3"

    [deps.Turing.extensions]
    TuringDynamicHMCExt = "DynamicHMC"
    TuringOptimExt = "Optim"

    [deps.Turing.weakdeps]
    DynamicHMC = "bbc10e6e-7c05-544b-b16e-64fede858acb"
    Optim = "429524aa-4258-5aef-a3af-852621145aeb"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.UnsafeAtomics]]
git-tree-sha1 = "6331ac3440856ea1988316b46045303bef658278"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.2.1"

[[deps.UnsafeAtomicsLLVM]]
deps = ["LLVM", "UnsafeAtomics"]
git-tree-sha1 = "323e3d0acf5e78a56dfae7bd8928c989b4f3083e"
uuid = "d80eeb9a-aca5-4d75-85e5-170c8b632249"
version = "0.1.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "7209df901e6ed7489fe9b7aa3e46fb788e15db85"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.65"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "5f24e158cf4cee437052371455fe361f526da062"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.6"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "801cbe47eae69adc50f36c3caec4758d2650741b"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.2+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "522b8414d40c4cbbab8dee346ac3a09f9768f25d"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.5+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "27798139afc0a2afa7b1824c206d5e87ea587a00"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.5"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "93284c28274d9e75218a416c65ec49d0e0fcdf3d"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.40+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╠═6975f6c3-8eec-47d9-ade5-88e8447b2dce
# ╠═14088e14-1fe6-41bd-a70a-38595fde2d41
# ╠═10cd1a4b-7190-4a1e-a867-80555d4bf730
# ╠═22f28b65-a270-4778-8933-57459e1249c7
# ╟─f3a514c8-bc3e-487f-af11-b43ff41c34e6
# ╟─e05e454f-c5b8-4f09-b745-a0265b5718a8
# ╠═a67ed97c-59ed-461b-803b-f2cb42498dda
# ╟─d6404e85-ceea-47c1-bc0e-c7c3deeb0ad6
# ╠═dcfef783-e51f-400b-9702-9068b8ef45b7
# ╟─62cf8524-0b46-434f-b813-f7befa9eb7f9
# ╠═079e0242-0d68-4e72-be44-9e317733071c
# ╠═e88ff706-532d-4890-ad6d-928aed3da558
# ╟─addadf8a-a757-4927-9e78-0b64a9bcf3eb
# ╠═247d18b1-dd98-49db-9185-2f76dae4b346
# ╟─c5e308ce-2a23-4773-b8d1-9dfb6a011d16
# ╠═35afca2b-8537-4ff1-a1c1-3021f3af3d4c
# ╟─38668cc3-5bcb-42fe-b89f-fed08dfacce5
# ╠═0195e6da-43f9-4ed3-bcad-7bee87fe23cf
# ╟─3172be0c-e3f6-43c3-b740-137584b19bcf
# ╠═1322553d-53ba-4166-a65f-a1350158dd46
# ╠═ae02befe-b6b1-49d4-8f85-00f59cd57f92
# ╟─140245a2-38cc-440a-96f2-7e375b60ef7b
# ╠═93cdbcce-1958-4476-a32a-9b1cf3ab9f17
# ╟─b5dcad20-acc6-446b-8231-cf9d6eb6bf6d
# ╠═1430a43d-d75d-4aff-8f43-a3c5e047c6eb
# ╠═a86c573b-d227-46f7-bd5b-81c44ef2b326
# ╠═e172123a-8186-4038-b9d6-edd5deffefd9
# ╟─776b5efd-f42b-49bc-98ed-4893eb216e2c
# ╠═9780cf8c-ca3f-486e-bae3-ca480dfd8d72
# ╟─b9d8d506-d88c-48a6-9e39-d33f0984047c
# ╠═ca8c5830-cd5d-421f-af9e-a0da2656d953
# ╠═a80e763b-9f3b-44aa-9abc-26362b4618bd
# ╟─b6a5f2a8-6e45-47ba-be97-961cfeb57699
# ╠═e0fb292f-8aba-4ca9-ba0c-1791cf063623
# ╟─c778178e-a6ab-4071-90d7-9cf543ef2606
# ╠═da6198fa-c7f2-4277-98f7-d6ba2e092815
# ╟─14bc80d0-3528-4cda-9311-02af2c322fc4
# ╠═5dd96185-34b6-460e-8f7f-e34d0c5eb469
# ╟─ea6f1ad4-781a-47a6-9cee-785da25f5967
# ╠═8adfcdf2-3da5-4d30-9301-b9c2e63d5490
# ╟─4e297c9d-e181-41ce-b1be-d5a6669bb903
# ╠═c8bd4ef6-7a94-46c7-be9a-6943043a6670
# ╟─e80cb9db-b93b-4b43-a55a-ee547fbc9d48
# ╠═fa2b5e49-933f-4fba-886c-0cf99f81eab5
# ╟─47a558cd-abbb-4dcb-8130-625d433d94fa
# ╠═82744b52-b300-4821-b095-19229835c5d4
# ╟─9a9ca097-db59-4e94-b2d5-47f4dbde48a3
# ╠═ea03ef5e-6332-41e7-a07b-55109cb97697
# ╠═7750a022-ef28-4104-84ac-ba43c03efaf5
# ╟─ca353ff4-fabb-4c72-8ba4-4294def7a751
# ╠═3f4c5f34-216d-4b8e-84a1-2c6c8b83716d
# ╠═f7940bb6-00f5-4ea3-87a2-cdaf719f6147
# ╠═9d3de9ec-3b74-41fc-86cd-4fca5ed0d3fe
# ╠═a734af10-a0ff-49bf-87ab-a39c12ba68a2
# ╠═6e5e457f-e829-4bf4-95ee-cc71f5f05818
# ╠═2b6cf6d6-5d97-46f5-8359-c069342afa82
# ╟─19a1560d-3524-484d-a4ca-e1f9678f3f4e
# ╠═1f67c25c-1ef2-4d07-aa04-e5e8e06b6092
# ╟─c70321cd-ef51-41bd-9744-aeec09bf2efb
# ╠═297964d6-d5ff-46a1-8b96-c6bba81f6be4
# ╟─5ce38815-7ef8-44fc-899e-e2158d6ca291
# ╠═cbd9fc1e-6035-4513-b113-e36dd36c56b9
# ╠═a952f503-d7d7-4b31-8a5a-510770d15bea
# ╟─65c3888c-490f-46dd-9a41-c9fa135c5f3c
# ╠═216bb14d-f167-4b18-87e2-ae635a413622
# ╠═a021c167-1334-4aaf-89b6-1b363ea90978
# ╟─087b4dea-68af-4022-a35b-b2e09e68ca59
# ╠═ae7e3091-3917-4bc7-80d6-2770c72f374f
# ╟─e2403387-7305-4c58-98f8-c9e60061a23f
# ╠═51f74f4e-3b24-49ef-8221-0c74de9a8fd4
# ╟─d2ef4d6c-6c77-4b88-8443-dd3dea5f4fa4
# ╠═bf0f80c3-9049-4d0f-82cb-c2561fcc2aae
# ╠═0290f649-11df-4d3e-81af-626d196316f3
# ╠═8b3dc9ea-fdac-4080-b284-4a8d3cf16f28
# ╠═f60af3f5-1048-49fb-a3a7-3f0b2b2884f9
# ╟─67f5cc12-e357-4104-a3c6-7d3e381c281c
# ╠═540dcd1c-a44a-4c8a-9048-3cb5c382c278
# ╟─515ddb6a-5917-43c6-8312-96e57753b478
# ╠═151de157-f083-476b-8fb4-4f4814f0a062
# ╟─22c6ef5a-6a54-4be6-a3ff-3a13382d5bcf
# ╠═ef44deb9-b1f0-40a0-82a7-b970ed7fdfc0
# ╠═3ad96e80-177a-4b74-8d2f-9e1efc60ed4c
# ╠═7bb690a7-04dc-4f0e-aabb-13190ec2cc6a
# ╠═b9e52a64-ed27-49f0-b5d4-8f1d5af96457
# ╠═ab439d9e-d427-4262-bacc-4e88759a04f4
# ╠═0e79c84d-c016-4ecd-bea1-b4a80f84e007
# ╟─cd0a1208-2721-4629-9784-1b1e137cd966
# ╠═f768fbac-9982-44b8-9518-883760201de7
# ╟─bb371f48-2aaf-4205-a3ae-832da7254c16
# ╠═989be969-0dd0-4ef9-9203-df59af4b7825
# ╠═46c1aa71-6177-4d9c-9fab-996cfe4187b9
# ╠═239721f2-2dcf-478c-9266-ff186389185f
# ╟─c0e44ae1-e1c9-481e-a748-e94d4fa4932b
# ╠═3bc11427-36c1-4a75-b8dd-36a61a7e8350
# ╠═9b4c1ba4-74c3-4ed0-a558-f326940335ff
# ╠═7785e952-0275-44fd-a91a-fe08ee29368c
# ╠═baa3ad4e-539a-4f2a-9daa-c0183472ca66
# ╠═cd0358bb-f42e-4cd5-9dee-d0104cf3e07b
# ╟─0cafb1cb-956f-4ba2-afb6-3546c8f2d990
# ╟─511ee325-6aab-4cb3-8e55-33f838007758
# ╠═766883c2-ff90-4ca3-bab4-278942bdf308
# ╟─94195529-1d37-42cd-a5ee-9cf777b1f1bb
# ╠═8b4d82a3-8a9c-462f-a3cb-e5a71f94ec58
# ╟─b5d03baa-7e2c-403b-8d65-983dad4723f3
# ╠═a49dc7d3-4743-4464-aabc-7a7fceb74bf6
# ╠═db5ad34c-65c4-432d-a6e2-c420d43d7faa
# ╠═8e16a5cb-1cd3-4ad5-8297-af28e503aac8
# ╠═6d083875-75e2-47b9-b58f-6cb6b5d95172
# ╠═53acdf8b-a783-4272-bd5e-147a5b0e1ae8
# ╠═957254e4-c9c0-4954-b6c4-24e6146151ff
# ╠═6c57ae79-8980-426a-bc27-d18b92e094e2
# ╠═6e4d7d10-1dae-40f6-9c1b-8ee0924c2f16
# ╠═8d45b2c3-8305-477d-9f88-9c4152d726df
# ╠═99cee684-134c-4def-bb8c-c7cb03cd06ab
# ╟─13c4376e-616f-44dc-8185-87a4f2def727
# ╠═a5b54f8c-e670-4402-b218-8164170db469
# ╟─d0af8d51-7fc8-41cd-9a51-20775f11e1e1
# ╠═7d462c37-f172-4624-8012-02a44bb19a4e
# ╠═9e910804-5539-424f-a769-d22302020888
# ╠═f603a8ff-ba98-4bc0-80c6-a89243f83a3d
# ╠═36a6c6f8-3189-4cbe-a1de-dc93d02ccd51
# ╠═0646f466-3247-495a-8a7c-742256f64b11
# ╟─199010f9-bce8-4318-9912-3ed575ef21fe
# ╠═448b3506-728b-4631-b6cd-62df0cd39f45
# ╟─88c22cbd-5f7b-4e50-9a0e-0ab22b857559
# ╠═2b7adad6-be9f-45f8-9385-e5d31ed69067
# ╠═fddece02-6e0e-466c-97aa-bb96fa1fe013
# ╠═aff8e056-ba12-4e13-b6bf-63da9181b8b7
# ╠═421f0b89-bf68-4546-8184-445b2e2de9fe
# ╠═88648d82-1c06-4b8c-bc91-58518b707905
# ╠═3ab8e039-da2e-4cca-9a1c-94d929a7aa78
# ╟─ba34af6e-85ee-4ff1-95c1-f627405262a2
# ╠═4f6e6857-f4c7-4c28-94e5-d79390127cad
# ╟─6cb0c0a9-de92-44f8-8b76-931b8d7a8058
# ╠═c87d89c9-f559-4210-a300-cd42ffdc0779
# ╟─85909da0-d2b9-437e-ae47-44940fcd74df
# ╠═ac56b077-a1c9-4e6d-9c78-64a613df5bb0
# ╠═8d0e632f-03ea-439c-9f6e-d57f489fc92c
# ╟─0c1ed9fa-80b0-449a-83be-f9e208f5e536
# ╠═545a7bd1-e56c-45e9-bdd6-efd3e14632fc
# ╟─d83b7902-85fd-4424-8ffe-98f9da71d0ae
# ╟─d4fa0864-c56d-4665-89a6-e9197e037c56
# ╠═65ed2cfe-4404-4756-a9d5-3d55cae1eefd
# ╟─9d721d0d-abe9-4f42-bc98-6245aeb2b2d0
# ╟─028d2911-1c4a-4de1-ae4d-83e8fb153f12
# ╟─a9c9cc35-05f6-49c8-a3df-2f303cf7885c
# ╠═bc2d4336-7a19-4179-b791-e0cdb63cee0a
# ╟─ac7531ca-02b5-44e9-a418-f8c399e668e8
# ╠═7690e5a3-0473-48fa-963c-7b55d00a4105
# ╟─37e13f5d-0d82-4c2d-a1c2-479843eb7886
# ╠═3e8dd9ee-867e-43e5-b74c-d4f586bfd41b
# ╠═ed01aaf1-1c1e-4441-9999-1a04916df98e
# ╟─13d724a6-6902-4ec1-86e8-3c6a86e3957d
# ╠═aff70993-542e-4be7-99cb-160059a84e06
# ╠═dae83135-14da-4b46-b974-782d4ec4e4d9
# ╠═4f2e5271-e301-4f44-8474-ec9b5c97150c
# ╠═2f08dbb9-f31a-4fc2-ab11-fdfb88e5ddb9
# ╠═cc4700af-9c54-4060-a9ce-8c426e25eeb1
# ╠═cf9ade41-251c-4320-a606-4e877b996632
# ╠═d6f7eff0-9142-4635-81f3-11ec8c2a500b
# ╠═06cf33fc-f0f8-4cc7-b811-d035fbf40b20
# ╠═d6fa7523-c369-454d-bce9-a12132fe01f1
# ╠═e2fbc671-9ea3-4ba3-b847-b71ed7a67933
# ╠═17745763-875a-45f3-9b1a-8c4c80d5060c
# ╠═eca35612-759b-43ef-a94b-69dc4856317f
# ╠═040be4fa-9158-45da-a9cd-ac419d59bb56
# ╠═4b34baa2-9528-4311-a54d-33f080f64103
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
