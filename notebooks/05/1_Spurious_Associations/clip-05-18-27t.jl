### A Pluto.jl notebook ###
# v0.12.3

using Markdown
using InteractiveUtils

# ╔═╡ 0831131a-077d-11eb-2412-fbd0e7cc3436
using DrWatson

# ╔═╡ 08314e98-077d-11eb-1d86-f9c9dd52e2c3
begin
	@quickactivate "SR2TuringPluto"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ 041b4048-077d-11eb-20b8-6f8b8ee72626
md"## Clip-05-18-27t.jl"

# ╔═╡ f2973922-077b-11eb-3833-9f97491e0ae2
begin
	df = CSV.read(sr_datadir("WaffleDivorce.csv"), DataFrame)
	df.D = zscore(df.Divorce)
	df.M = zscore(df.Marriage)
	z_age, unz_age = zscore_transform(df.MedianAgeMarriage)
	df.A = z_age(df.MedianAgeMarriage)
end;

# ╔═╡ 3d995aa2-08e1-11eb-09f7-35b6aa59fe73
@model function m5_3_A(A, M, D)
    # A -> M
    aM ~ Normal(0, 0.2)
    bAM ~ Normal(0, 0.5)
    σ_M ~ Exponential(1)
    μ_M = lin(aM, A, bAM)
    M ~ MvNormal(μ_M, σ_M)
    # A -> D <- M
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    bA ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bA, M, bM)
    D ~ MvNormal(μ, σ)
end

# ╔═╡ 3d99976a-08e1-11eb-344e-f193f1387780
begin
	quap5_3_At = quap(m5_3_A(df.A, df.M, df.D))

	A_seq = range(-2, 2, length = 30)
	post5_3_At = DataFrame(rand(quap5_3_At.distr, 1000)', quap5_3_At.params)

	μM = lin(post5_3_At.aM', A_seq, post5_3_At.bAM')
	M = rand.(Normal.(μM, post5_3_At.σ_M'))
	μD = lin(post5_3_At.a', A_seq, post5_3_At.bA', M, post5_3_At.bM')
	D = rand.(Normal.(μD, post5_3_At.σ')) |> meanlowerupper

	plot(A_seq, D.mean, ribbon = (D.mean .- D.lower, D.upper .- D.mean))
	plot!(legend = false, xlabel = "manipulated A", ylabel = "counterfactual D")
end

# ╔═╡ 3dac58dc-08e1-11eb-33df-43cc24f8f0a5
begin
	M_seq = range(-2, 2, length = 30)
	A_seq1 = zeros(length(M_seq))
	M2 = rand.(Normal.(lin(post5_3_At.aM', A_seq, post5_3_At.bAM'), post5_3_At.σ_M'))
	μD2 = lin(post5_3_At.a', A_seq1, post5_3_At.bA', M_seq, post5_3_At.bM')
	D1 = rand.(Normal.(μD2, post5_3_At.σ')) |> meanlowerupper

	plot(M_seq, D1.mean, ribbon = (D1.mean .- D1.lower, D1.upper .- D1.mean))
	plot!(legend = false, xlabel = "manipulated M", ylabel = "counterfactual D")
end

# ╔═╡ 3dad5228-08e1-11eb-1b14-6d412246f714
begin
	post = DataFrame(rand(quap5_3_At.distr, 1000)', quap5_3_At.params)
	M_sim = rand.(Normal.(post.aM' .+ A_seq .* post.bAM', post.σ_M'))
	D_sim = rand.(Normal.(post.a' .+ A_seq .* post.bA', post.σ'))
end

# ╔═╡ 3d9a2388-08e1-11eb-28d1-131c48d85129
begin
	M1 = rand.(Normal.(lin(post.aM', (20.0, 30.0) |> z_age, post5_3_At.bAM'), post5_3_At.σ_M'))
	μD1 = lin(post5_3_At.a', (20.0, 30.0) |> z_age, post5_3_At.bA', M1, post5_3_At.bM') |> meanlowerupper
	μD1.mean[2] - μD1.mean[1]
end

# ╔═╡ 3db8aeb6-08e1-11eb-3947-8fb3cd806df2
md"## End of clip-05--19-27t.jl"

# ╔═╡ Cell order:
# ╟─041b4048-077d-11eb-20b8-6f8b8ee72626
# ╠═0831131a-077d-11eb-2412-fbd0e7cc3436
# ╠═08314e98-077d-11eb-1d86-f9c9dd52e2c3
# ╠═f2973922-077b-11eb-3833-9f97491e0ae2
# ╠═3d995aa2-08e1-11eb-09f7-35b6aa59fe73
# ╠═3d99976a-08e1-11eb-344e-f193f1387780
# ╠═3d9a2388-08e1-11eb-28d1-131c48d85129
# ╠═3dac58dc-08e1-11eb-33df-43cc24f8f0a5
# ╠═3dad5228-08e1-11eb-1b14-6d412246f714
# ╟─3db8aeb6-08e1-11eb-3947-8fb3cd806df2
