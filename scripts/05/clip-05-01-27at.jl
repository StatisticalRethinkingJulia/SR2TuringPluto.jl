
using Markdown
using InteractiveUtils

using DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Clip-05-01-27t.jl"

begin
	df = CSV.read(sr_datadir("WaffleDivorce.csv"), DataFrame)
	df.D = zscore(df.Divorce)
	df.M = zscore(df.Marriage)
	df.A = zscore(df.MedianAgeMarriage)
end

std(df.MedianAgeMarriage)

@model function m5_1_A(A, D)
    a ~ Normal(0, 0.2)
    bA ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bA)
    D ~ MvNormal(μ, σ)
end

begin
	m5_1_At = m5_1_A(df.A, df.D)
	prior5_1_At = sample(m5_1_At, Prior(), 50) |> DataFrame
	prior5_1_At = prior5_1_At[:, 3:end]
	Text(precis(prior5_1_At; io=String))
end

begin
	x = -2:0.1:2
	plot()
	for r in eachrow(prior5_1_At)
		p = lin(r.a, x, r.bA)
		plot!(x, p, color = :black, alpha = 0.4)
	end
	plot!(legend = false)

	quap5_1_At = quap(m5_1_At)
	dfa5_1_At = DataFrame(rand(quap5_1_At.distr, 1000)', quap5_1_At.params)

	A_seq = range(-3, 3.2, length = 30)
	mu5_1_At = lin(dfa5_1_At.a', A_seq, dfa5_1_At.bA') |> meanlowerupper

	scatter(df.A, df.D, alpha = 0.4, legend = false)
	plot!(A_seq, mu5_1_At.mean, ribbon =
		(mu5_1_At.mean .- mu5_1_At.lower, mu5_1_At.upper .- mu5_1_At.mean))
	vline!([0])
end

@model function m5_1_M(M, D)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, M, bM)
    D ~ MvNormal(μ, σ)
end

quap5_1_Mt = quap(m5_1_M(df.M, df.D))

begin
	dfa5_1_Mt = DataFrame(rand(quap5_1_Mt.distr, 1000)', quap5_1_Mt.params)
	M_seq = range(-3, 3.2, length = 30)
	mu5_1_Mt = lin(dfa5_1_Mt.a', M_seq, dfa5_1_Mt.bM') |> meanlowerupper
	scatter(df.M, df.D, alpha = 0.4, legend = false)
	plot!(M_seq, mu5_1_Mt.mean, ribbon = 
		(mu5_1_Mt.mean .- mu5_1_Mt.lower, mu5_1_Mt.upper .- mu5_1_Mt.mean))
	vline!([0])
end

@model function m5_1_A_M(A, M, D)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    bA ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bA, M, bM)
    D ~ MvNormal(μ, σ)
end

begin
	m5_3_A_Mt = m5_1_A_M(df.A, df.M, df.D)
	quap5_1_A_Mt = quap(m5_3_A_Mt)
end

xm = [
    quap5_1_At.coef.bA,       # bA
    NaN,
    quap5_1_A_Mt.coef.bA,
    NaN,                	  # Spaceholder
    NaN,                	  # bM
    quap5_1_Mt.coef.bM,
    quap5_1_A_Mt.coef.bM,
];

xerr = sqrt.([
    quap5_1_At.vcov[2, 2],    # bA
    NaN,
    quap5_1_A_Mt.vcov[3, 3],
    NaN,                      # Spaceholder
    NaN,                      # bM
    quap5_1_Mt.vcov[2, 2],
    quap5_1_A_Mt.vcov[2, 2],
]);

ylab = [
    "bA m5.1_At",
    "bA m5.1_Mt",
    "bA m5.1_A_Mt",
    "",
    "bM m5.1_At",
    "bM m5.1_Mt",
    "bM m5.1_A_Mt"
];

scatter(x, ylab, xerr = xerr, legend = false)



N = 50
age = randn(N)
mar = rand.(Normal.(-age))
div = rand.(Normal.(age))

@model function divorce_AM(A, M)
    a ~ Normal(0, 0.2)
    bAM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bAM)
    M ~ MvNormal(μ, σ)
end

sort!(d, :A)  # for nice looking graphs

q5_4 = quap(divorce_AM(d.A, d.M))
post = DataFrame(rand(q5_4.distr, 1000)', q5_4.params)
mu = lin(post.a', d.A, post.bAM') |> meanlowerupper
resid = d.M .- mu.mean

scatter(d.A, d.M, legend = false)
plot!(d.A, mu.mean, ribbon = (mu.mean .- mu.lower, mu.upper .- mu.mean))

post = DataFrame(rand(q5_3.distr, 1000)', q5_3.params)
μ = lin(post.a', d.A, post.bA', d.M, post.bM') |> meanlowerupper

scatter(d.D, μ.mean, yerror = (μ.mean .- μ.lower, μ.upper .- μ.mean), legend = false)
plot!(identity, range(extrema(d.D)..., length = 10))
plot!(xlabel = "Observed divorce", ylabel = "Predicted divorce")

xy = [d.D μ.mean]
for loc in ("ID", "UT", "ME", "RI")
    coord = xy[d.Loc .== loc, :]
    annotate!(coord..., loc, :bottom)
end
plot!()

N = 100
x_real = randn(N)
x_spur = rand.(Normal.(x_real))
y = rand.(Normal.(x_real))
d = DataFrame((; y, x_real, x_spur))  # the semicolon turns it into a NamedTuple which then
                                      # gives the DataFrame the names of the columns

@df d corrplot([:y :x_real :x_spur])  # either
corrplot(Matrix(d), label = names(d)) # or

d = DataFrame(CSV.File(datadir("exp_raw/WaffleDivorce.csv")))
d.D = zscore(d.Divorce)
d.M = zscore(d.Marriage)
z_age, unz_age = zscore_transform(d.MedianAgeMarriage)
d.A = z_age(d.MedianAgeMarriage)

@model function divorce_19(A, M, D)
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

q5_3_A = quap(divorce_19(d.A, d.M, d.D))

A_seq = range(-2, 2, length = 30)
post = DataFrame(rand(q5_3_A.distr, 1000)', q5_3_A.params)

μM = lin(post.aM', A_seq, post.bAM')
M = rand.(Normal.(μM, post.σ_M'))
μD = lin(post.a', A_seq, post.bA', M, post.bM')
D = rand.(Normal.(μD, post.σ')) |> meanlowerupper

plot(A_seq, D.mean, ribbon = (D.mean .- D.lower, D.upper .- D.mean))
plot!(legend = false, xlabel = "manipulated A", ylabel = "counterfactual D")

M = rand.(Normal.(lin(post.aM', (20.0, 30.0) |> z_age, post.bAM'), post.σ_M'))
μD = lin(post.a', (20.0, 30.0) |> z_age, post.bA', M, post.bM') |> meanlowerupper
μD.mean[2] - μD.mean[1]

M_seq = range(-2, 2, length = 30)
A_seq = zeros(length(M_seq))
M = rand.(Normal.(lin(post.aM', A_seq, post.bAM'), post.σ_M'))
μD = lin(post.a', A_seq, post.bA', M_seq, post.bM')
D = rand.(Normal.(μD, post.σ')) |> meanlowerupper

plot(M_seq, D.mean, ribbon = (D.mean .- D.lower, D.upper .- D.mean))
plot!(legend = false, xlabel = "manipulated M", ylabel = "counterfactual D")

A_seq = range(-2, 2, length = 30)
post = DataFrame(rand(q5_3_A.distr, 1000)', q5_3_A.params)
M_sim = rand.(Normal.(post.aM' .+ A_seq .* post.bAM', post.σ_M'))
D_sim = rand.(Normal.(post.a' .+ A_seq .* post.bA', post.σ'))



