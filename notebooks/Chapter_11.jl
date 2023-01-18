### A Pluto.jl notebook ###
# v0.19.3

using Markdown
using InteractiveUtils

# ╔═╡ 14088e14-1fe6-41bd-a70a-38595fde2d41
using Pkg, DrWatson

# ╔═╡ 10cd1a4b-7190-4a1e-a867-80555d4bf730
begin
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
	PRECIS(prior11_2t_df)
end

# ╔═╡ 93cdbcce-1958-4476-a32a-9b1cf3ab9f17
let
	f = i -> @. logistic(prior11_2t_df.a + prior11_2t_df[!,"b[$i]"])
	p1 = map(f, 1:4)
	density(abs.(p1[1] .- p1[2]), bandwidth=0.01)
end

# ╔═╡ 140245a2-38cc-440a-96f2-7e375b60ef7b
md" #### Code 11.8"

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
	PRECIS(m11_4_df)
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
	PRECIS(m11_5_df)
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
	PRECIS(m11_7_df)
end

# ╔═╡ 67f5cc12-e357-4104-a3c6-7d3e381c281c
md" #### Code 11.30"

# ╔═╡ 540dcd1c-a44a-4c8a-9048-3cb5c382c278
let
	diff_a = m11_7_df."a[1]" .- m11_7_df."a[2]"
	diff_p = @. logistic(m11_7_df."a[1]") - logistic(m11_7_df."a[2]")
	PRECIS(DataFrame(diff_a=diff_a, diff_p=diff_p))
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
	PRECIS(m11_8_df)
end

# ╔═╡ ab439d9e-d427-4262-bacc-4e88759a04f4
md" #### Code 11.33"

# ╔═╡ 0e79c84d-c016-4ecd-bea1-b4a80f84e007
let
	diff_a = m11_8_df."a[1]" .- m11_8_df."a[2]"
	diff_p = logistic.(m11_8_df."a[1]") .- logistic.(m11_8_df."a[2]")
	PRECIS(DataFrame(diff_a=diff_a, diff_p=diff_p))
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
	PRECIS(m11_8_df)
end

# ╔═╡ c0e44ae1-e1c9-481e-a748-e94d4fa4932b
md" #### Code 11.33"

# ╔═╡ 3bc11427-36c1-4a75-b8dd-36a61a7e8350
begin
	diff_a = m11_8_df."a[1]" .- m11_8_df."a[2]"
	diff_p = logistic.(m11_8_df."a[1]") .- logistic.(m11_8_df."a[2]")
	PRECIS(DataFrame(diff_a=diff_a, diff_p=diff_p))
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
PRECIS(m11_11_df)

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
	PRECIS(m11_12_df)
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
	PRECIS(DataFrame(λ_old=λ_old, λ_new=λ_new))
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
	PRECIS(m11_13_df)
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
	PRECIS(DataFrame(p_diff=p_diff))
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
	PRECIS(m11_14_df)
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
