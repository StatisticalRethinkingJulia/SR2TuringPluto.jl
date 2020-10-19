
using Markdown
using InteractiveUtils

using DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Clip-05-18-27t.jl"

begin
	df = CSV.read(sr_datadir("WaffleDivorce.csv"), DataFrame)
	df.D = zscore(df.Divorce)
	df.M = zscore(df.Marriage)
	z_age, unz_age = zscore_transform(df.MedianAgeMarriage)
	df.A = z_age(df.MedianAgeMarriage)
end;

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

begin
	M_seq = range(-2, 2, length = 30)
	A_seq1 = zeros(length(M_seq))
	M2 = rand.(Normal.(lin(post5_3_At.aM', A_seq, post5_3_At.bAM'), post5_3_At.σ_M'))
	μD2 = lin(post5_3_At.a', A_seq1, post5_3_At.bA', M_seq, post5_3_At.bM')
	D1 = rand.(Normal.(μD2, post5_3_At.σ')) |> meanlowerupper

	plot(M_seq, D1.mean, ribbon = (D1.mean .- D1.lower, D1.upper .- D1.mean))
	plot!(legend = false, xlabel = "manipulated M", ylabel = "counterfactual D")
end

begin
	post = DataFrame(rand(quap5_3_At.distr, 1000)', quap5_3_At.params)
	M_sim = rand.(Normal.(post.aM' .+ A_seq .* post.bAM', post.σ_M'))
	D_sim = rand.(Normal.(post.a' .+ A_seq .* post.bA', post.σ'))
end

begin
	M1 = rand.(Normal.(lin(post.aM', (20.0, 30.0) |> z_age, post5_3_At.bAM'), post5_3_At.σ_M'))
	μD1 = lin(post5_3_At.a', (20.0, 30.0) |> z_age, post5_3_At.bA', M1, post5_3_At.bM') |> meanlowerupper
	μD1.mean[2] - μD1.mean[1]
end

md"## End of clip-05--19-27t.jl"

