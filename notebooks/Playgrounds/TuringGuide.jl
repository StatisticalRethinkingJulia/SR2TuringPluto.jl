### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 2d3d462e-d6bf-4ebd-ab34-db38f4261a11
using Pkg, DrWatson

# ╔═╡ 18b5c138-9324-4c5a-86ca-be1e825ed11e
begin
	using PlutoUI
	using Random
	using Distributions
	using StatsPlots
	using StatsBase
	using LaTeXStrings
	using CSV
	using DataFrames
	using PrettyTables
	using StatisticalRethinking
	using StatisticalRethinkingPlots
	using LinearAlgebra
	using Logging
	using AxisKeys
	using NamedDims
	using NamedArrays
	using DataStructures
	using Random
	using Optim
	using Turing
	using ParetoSmooth
end

# ╔═╡ 58b9da10-8cbd-42e4-bddf-ec4995ac756c
md" ## Turing guide"

# ╔═╡ d3c6c774-a158-485c-acd1-e366c24df4bc
md" ##### Widen notebook."

# ╔═╡ 0f75e85c-4506-483f-b0f8-a31ee698bf04
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

# ╔═╡ cb180e21-3778-40f1-93d8-af442d1be559
md" ##### Commonly used packages in SR2TuringPluto"

# ╔═╡ 1bc33245-85cf-4ec9-be87-9cde44527cbe
md"

!!! note

In most cases I like to enable below line. But if ParetoSmooth is needed, it will hold back many newer package versions. Thus ParetoSmooth.jl is not included in SR2TuringPluto.jl.

"

# ╔═╡ 2d6ec605-7efa-48a9-8e69-4e83294c716e
#Pkg.activate(expanduser("~/.julia/dev/SR2TuringPluto"))

# ╔═╡ a2d950f8-d780-4abe-af37-b189bd4a5e88
md"##### Setting default attributes for logging and plots."

# ╔═╡ bbb4a889-77e2-449b-9258-883e7434e641
begin
	Logging.disable_logging(Logging.Warn);
end;

# ╔═╡ 496b257c-ea62-40b2-b1c9-e92fdda1d688
begin
	howell = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	howell1 = howell[howell.age .>= 18, :]
	describe(howell)
end

# ╔═╡ 96760a8e-8561-4576-a7c7-70596fa20d33
md" ### Statistical Rethinking (SR2) models."

# ╔═╡ b83bc4d0-8c65-40a8-9187-70ef2a0ff38d
md"In the book and associated R package `rethinking`, statistical models are defined as illustrated below:

```
flist <- alist(
  height ~ dnorm( mu , sigma ) ,
  mu <- a + b*weight ,
  a ~ dnorm( 156 , 100 ) ,
  b ~ dnorm( 0 , 10 ) ,
  sigma ~ dunif( 0 , 50 )
)
```"

# ╔═╡ b583eb23-4045-420c-80d5-4f339ebb6ad4
md" ##### An equivalent Turing model:"

# ╔═╡ d5d66e16-487e-4214-b702-728dc3db32cc
@model function ppl4_1(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* weights
    for i in eachindex(heights)
        heights[i] ~ Normal(μ[i], σ)
    end
end

# ╔═╡ b14217e8-1f70-44cc-8d3c-528ab55e27ef
md"
!!! note

The sequence of the statements matter in Turing models!"

# ╔═╡ b3adef33-c111-4d20-b277-a5346eae9f23
begin
	m4_1t = ppl4_1(howell1.weight, howell1.height)
	post4_1t = sample(m4_1t, NUTS(), 1000)
	post4_1t_df = DataFrame(post4_1t)
	post4_1t_df
end

# ╔═╡ 6869d583-1590-4210-b696-b63fd7857e09
plot(post4_1t)

# ╔═╡ 1f4e22a4-5fba-4331-9bf9-dd886ce45e74
begin
	x = 30.0:0.1:64.0
	preds = mean(post4_1t[:a]) .+ mean(post4_1t[:b]) .* x
	plot(x, preds, color=:darkblue, label="Regression line")
	scatter!(howell1.weight, howell1.height, marker=:cross, markersize=3,
		color=:darkred, label="Observations")
end

# ╔═╡ cbe1028d-829d-4ed4-840d-d74f151f55bb
CHNS(post4_1t[[:a, :b, :σ, :lp, :log_density]])

# ╔═╡ 75ee3adf-2988-4f42-a0b5-21abe5829120
md" ### Priors of Turing models."

# ╔═╡ 63db6c9f-ba40-41a2-975d-798d9f2cbeab
let
	priors4_1t = sample(m4_1t, Prior(), 1000)
	DataFrame(priors4_1t)
end

# ╔═╡ 534edcc4-2dac-4ddb-a6a9-55ee38a3e7ae
md" ### Quadratic approximation."

# ╔═╡ 7c8e8a61-b1c1-4f9a-addc-61a987c2c5d2
map_estimate = optimize(m4_1t, MAP())

# ╔═╡ a9d5118d-8fc7-4353-bb20-6b87ae0e9247
begin
	â, b̂, σ̂ = map_estimate.values
	(â = â, b̂ = b̂, σ̂ = σ̂)
end

# ╔═╡ f41d6486-7589-4b67-9eaf-92956622226a
let
	x = 30:0.1:65
	plot(x, â .+ b̂ .* x)
	scatter!(howell1.weight, howell1.height, leg=false)
end

# ╔═╡ 9f92d0a7-45ce-4148-b6a5-7ede9dc09234
md" ### Prediction."

# ╔═╡ 1cc80dac-1120-4520-b73f-0e31dc1d8721
begin
	x_test = [31, 40, 50, 60, 63]
	m_test = ppl4_1(x_test, fill(missing, length(x_test)))
	pred4_1t = predict(m_test, post4_1t)
	describe(DataFrame(pred4_1t))
end

# ╔═╡ 159e07b6-dcd7-4ca7-9d3e-f5059b1b8803
md"##### Get the mean predicted values."

# ╔═╡ 85121bc0-447a-4b79-ac6c-ce91871ca7d8
begin
	pred_summary_df = DataFrame()
	pred_summary_df.parameters = names(pred4_1t)
	pred_summary_df.mean = mean(pred4_1t)[:, 2]
	pred_summary_df.std = vec(std(Array(group(pred4_1t, :heights)); dims = 1))
	pred_summary_df.lower = hpd(pred4_1t; alpha=0.11)[:, 2]
	pred_summary_df.upper = hpd(pred4_1t; alpha=0.11)[:, 3]
	pred_summary_df
end

# ╔═╡ 84b6ab83-7e65-477f-b87f-298d8b5eacea
let
	scatter(howell1.weight, howell1.height, lab="Observations",leg=:bottomright)
	scatter!(x_test, pred_summary_df.mean,
    	markersize=5, color=:red, lab="Predictions at x_test")
	x = range(minimum(howell1.weight), stop=maximum(howell1.weight), length=200)
	g(x) = mean(post4_1t[:a]) + mean(post4_1t[:b]) * x
	plot!(x, g.(x), lab="Regression line", leg=:bottomright)
end

# ╔═╡ 2d373995-f303-48b7-8cc3-1daf35d5ed36
md" ### StatisticalRethinking link() function."

# ╔═╡ 91a627e0-a46c-40f8-a5c2-2f3b25187074
@bind N Slider(1:1000, default=500)

# ╔═╡ a492bfeb-452c-4138-b5a2-6d2c4cdd114b
let
	xseq = 31:63
	μ = StatisticalRethinking.link(post4_1t_df, [:a, :b], xseq)
	μ = hcat(μ...)

	p = plot(; xlim=xseq, ylim=xseq, 
		xlab="weight", ylab="height", 
		title="Regression lines ($N) from draws",
		leg=false
	)
	for y ∈ first(eachrow(μ), N)
		plot!(p, xseq, y; c=:lightgrey, alpha=0.3)
	end
	p
end

# ╔═╡ ccf333ff-af94-439a-b41c-5be04c076007
md"""
!!! note

Click on "Live docs" and place cursor on link to see more help. 

Click little down arrow to the right to remove live docs again.
"""

# ╔═╡ 65eaad30-7f58-4eb1-b458-07688dfeac81
md" ### Model comparison"

# ╔═╡ 3283d295-c60e-44bd-9547-8285cdbce9a9
Random.seed!(129111)

# ╔═╡ de4e3b12-174b-4740-acab-db6d2a55eb0a
begin
	df = CSV.read(sr_datadir("WaffleDivorce.csv"), DataFrame)
	df.D = zscore(df.Divorce)
	df.M = zscore(df.Marriage)
	df.A = zscore(df.MedianAgeMarriage)
	data = (D=df.D, A=df.A)
end;

# ╔═╡ 28bf03e4-2899-4f34-a8aa-892937ede3f9
function lin(a, b, c, x...)
	result = @. a + b * c
	for i in 1:2:length(x)
		@. result += x[i] * x[i + 1]
	end
	return result
end

# ╔═╡ 30707462-75f5-4da9-9677-d7208f20a4fb
@model function ppl5_1(A, D)
	a ~ Normal(0, 0.2)
	bA ~ Normal(0, 0.5)
	σ ~ Exponential(1)
	for i in eachindex(D)
		μ = lin(a, A[i], bA)
		D[i] ~ Normal(μ, σ)
	end
end

# ╔═╡ 18b3a4a5-2506-4ced-8917-3da988c94eba
begin
	m5_1t = ppl5_1(df.A, df.D)
	chn5_1t = sample(m5_1t, NUTS(1000, 0.9), 1000)
	describe(chn5_1t)
end

# ╔═╡ 91b50108-414d-4520-8e57-5eaa58061b1a
@model function ppl5_2(M, D)
	a ~ Normal(0, 0.2)
	bM ~ Normal(0, 0.5)
	σ ~ Exponential(1)
	for i in eachindex(D)
		μ = lin(a, M[i], bM)
		D[i] ~ Normal(μ, σ)
	end
end

# ╔═╡ f11d0597-8962-4cd7-8637-aef99abce80a
begin
	m5_2t = ppl5_2(df.M, df.D)
	chn5_2t = sample(m5_2t, NUTS(1000, 0.9), 1000)
	DataFrame(chn5_2t)
end

# ╔═╡ a0834a65-1f5e-497a-89ab-36adf38d1118
@model function ppl5_3(A, M, D)
	a ~ Normal(0, 0.2)
	bA ~ Normal(0, 0.5)
	bM ~ Normal(0, 0.5)
	σ ~ Exponential(1)
	for i in eachindex(D)
		μ = a + M[i] * bM + A[i] * bA
		D[i] ~ Normal(μ, σ)
	end
end

# ╔═╡ ba0b9017-12c9-446d-a8dd-647d8735d2b0
begin
	m5_3t = ppl5_3(df.A, df.M, df.D)
	chn5_3t = sample(m5_3t, NUTS(1000, 0.9), 1000)
	DataFrame(chn5_3t)
end

# ╔═╡ e0febcf9-6e87-4ecc-8ef9-e27ab79af2fc
function model_elpd(psis::T; digits=2) where {T <: PsisLoo}
	NamedArray(
		round.(Array(psis.estimates); digits=digits), 
		(OrderedDict(:cv_elpd=>1, :naive_lpd=>2, :p_eff=>3), 
		OrderedDict(:total=>1, :se_total=>2, :mean=>3, :se_mean=>4)),
               ("Statistic", "Value")
	)
end

# ╔═╡ c82d6f3b-6a27-4cff-9043-c301aa09d53e
begin
	#pw_lls5_1t = pointwise_log_likelihoods(m5_1t, chn5_1t)
	psis5_1t = psis_loo(m5_1t, chn5_1t)
	model_elpd(psis5_1t; digits=1)
end

# ╔═╡ 74977e3a-b168-4bcd-a260-7191c60130ac
begin
	#pw_lls5_2t = pointwise_log_likelihoods(m5_2t, chn5_2t)
	psis5_2t = psis_loo(m5_2t, chn5_2t)
	model_elpd(psis5_2t; digits=1)
end

# ╔═╡ 7d6161c6-4805-4ee9-8d9a-d8f95faf47d7
begin
	#pw_lls5_3t = pointwise_log_likelihoods(m5_3t, chn5_3t)
	psis5_3t = psis_loo(m5_3t, chn5_3t)
	model_elpd(psis5_3t; digits=1)
end

# ╔═╡ d8eb2079-2915-46e5-ac55-6849143aa194
md"""

```
df_waic = compare([m5_1s, m5_2s, m5_3s], :waic)
3×8 DataFrame
 Row │ models  WAIC     lppd     SE       dWAIC    dSE      pWAIC    weight  
     │ String  Float64  Float64  Float64  Float64  Float64  Float64  Float64 
─────┼───────────────────────────────────────────────────────────────────────
   1 │ m5.1s     125.6   118.49    12.56      0.0     0.0      3.57     0.73
   2 │ m5.3s     127.6   118.13    12.63      2.0     0.74     4.72     0.27
   3 │ m5.2s     139.2   133.4      9.8      13.6     9.15     2.88     0.0

compare([m5_1s, m5_2s, m5_3s], :psis)
3×8 DataFrame
 Row │ models  PSIS     lppd     SE       dPSIS    dSE      pPSIS    weight  
     │ String  Float64  Float64  Float64  Float64  Float64  Float64  Float64 
─────┼───────────────────────────────────────────────────────────────────────
   1 │ m5.1s     125.1   118.49    12.35      0.0     0.0      3.57     0.71
   2 │ m5.3s     126.9   118.13    12.42      1.8     0.72     4.72     0.29
   3 │ m5.2s     138.7   133.4      9.7      13.6     8.87     2.88     0.0
```
"""

# ╔═╡ 5e18ac09-76e1-44a1-93a9-408d6971a5f6
function compare_models(nt; digits=2)
	comps = loo_compare(nt)
	estimates = comps.estimates
	models = Pair{Symbol, Int}[]
	for (indx, name) in enumerate(estimates.model)
		append!(models, [name => indx])
	end
	
	NamedArray(
		round.(Array(estimates); digits=2), 
		(OrderedDict(models...), 
		OrderedDict(:cv_elpd=>1, :cv_avg=>2, :weight=>3)),
               ("Statistic", "Comparison")
	)
end

# ╔═╡ cc0b2fad-d508-4135-8203-72e3d130fea1
mc = compare_models((m5_1t=psis5_1t, m5_2t=psis5_2t, m5_3t=psis5_3t))

# ╔═╡ 3d2862bc-c217-4e8d-bf88-24e7e3b3eaa7
mc[:m5_3t, :weight]

# ╔═╡ 222a9455-f017-47ff-bd52-42179f4a0320
md" ### Graphs and DAGs."

# ╔═╡ faeea7c5-d921-4dbe-9f72-d3e0ae17748f
md" ### Working with ... ."

# ╔═╡ Cell order:
# ╟─58b9da10-8cbd-42e4-bddf-ec4995ac756c
# ╟─d3c6c774-a158-485c-acd1-e366c24df4bc
# ╠═0f75e85c-4506-483f-b0f8-a31ee698bf04
# ╟─cb180e21-3778-40f1-93d8-af442d1be559
# ╠═2d3d462e-d6bf-4ebd-ab34-db38f4261a11
# ╟─1bc33245-85cf-4ec9-be87-9cde44527cbe
# ╠═2d6ec605-7efa-48a9-8e69-4e83294c716e
# ╠═18b5c138-9324-4c5a-86ca-be1e825ed11e
# ╟─a2d950f8-d780-4abe-af37-b189bd4a5e88
# ╠═bbb4a889-77e2-449b-9258-883e7434e641
# ╠═496b257c-ea62-40b2-b1c9-e92fdda1d688
# ╟─96760a8e-8561-4576-a7c7-70596fa20d33
# ╟─b83bc4d0-8c65-40a8-9187-70ef2a0ff38d
# ╟─b583eb23-4045-420c-80d5-4f339ebb6ad4
# ╠═d5d66e16-487e-4214-b702-728dc3db32cc
# ╟─b14217e8-1f70-44cc-8d3c-528ab55e27ef
# ╠═b3adef33-c111-4d20-b277-a5346eae9f23
# ╠═6869d583-1590-4210-b696-b63fd7857e09
# ╠═1f4e22a4-5fba-4331-9bf9-dd886ce45e74
# ╠═cbe1028d-829d-4ed4-840d-d74f151f55bb
# ╟─75ee3adf-2988-4f42-a0b5-21abe5829120
# ╠═63db6c9f-ba40-41a2-975d-798d9f2cbeab
# ╟─534edcc4-2dac-4ddb-a6a9-55ee38a3e7ae
# ╠═7c8e8a61-b1c1-4f9a-addc-61a987c2c5d2
# ╠═a9d5118d-8fc7-4353-bb20-6b87ae0e9247
# ╠═f41d6486-7589-4b67-9eaf-92956622226a
# ╟─9f92d0a7-45ce-4148-b6a5-7ede9dc09234
# ╠═1cc80dac-1120-4520-b73f-0e31dc1d8721
# ╟─159e07b6-dcd7-4ca7-9d3e-f5059b1b8803
# ╠═85121bc0-447a-4b79-ac6c-ce91871ca7d8
# ╠═84b6ab83-7e65-477f-b87f-298d8b5eacea
# ╟─2d373995-f303-48b7-8cc3-1daf35d5ed36
# ╠═91a627e0-a46c-40f8-a5c2-2f3b25187074
# ╠═a492bfeb-452c-4138-b5a2-6d2c4cdd114b
# ╟─ccf333ff-af94-439a-b41c-5be04c076007
# ╟─65eaad30-7f58-4eb1-b458-07688dfeac81
# ╠═3283d295-c60e-44bd-9547-8285cdbce9a9
# ╠═de4e3b12-174b-4740-acab-db6d2a55eb0a
# ╠═28bf03e4-2899-4f34-a8aa-892937ede3f9
# ╠═30707462-75f5-4da9-9677-d7208f20a4fb
# ╠═18b3a4a5-2506-4ced-8917-3da988c94eba
# ╠═91b50108-414d-4520-8e57-5eaa58061b1a
# ╠═f11d0597-8962-4cd7-8637-aef99abce80a
# ╠═a0834a65-1f5e-497a-89ab-36adf38d1118
# ╠═ba0b9017-12c9-446d-a8dd-647d8735d2b0
# ╠═e0febcf9-6e87-4ecc-8ef9-e27ab79af2fc
# ╠═c82d6f3b-6a27-4cff-9043-c301aa09d53e
# ╠═74977e3a-b168-4bcd-a260-7191c60130ac
# ╠═7d6161c6-4805-4ee9-8d9a-d8f95faf47d7
# ╟─d8eb2079-2915-46e5-ac55-6849143aa194
# ╠═5e18ac09-76e1-44a1-93a9-408d6971a5f6
# ╠═cc0b2fad-d508-4135-8203-72e3d130fea1
# ╠═3d2862bc-c217-4e8d-bf88-24e7e3b3eaa7
# ╟─222a9455-f017-47ff-bd52-42179f4a0320
# ╟─faeea7c5-d921-4dbe-9f72-d3e0ae17748f
