### A Pluto.jl notebook ###
# v0.11.12

using Markdown
using InteractiveUtils

# ╔═╡ a2e13c44-f08e-11ea-13df-53d884942d7e
using Pkg, DrWatson

# ╔═╡ 021cabc4-f08f-11ea-1fc2-f19729723f2c
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StanSample
	using StatisticalRethinking
end

# ╔═╡ ca40913e-f093-11ea-28be-15a831f7c1c2
md"## Clip-02-06-07t.jl"

# ╔═╡ 0c39592c-f094-11ea-3cb9-17977c11b058
md"##### Pull in Stansample to compare the quap results"

# ╔═╡ 7554b2f8-f08f-11ea-38c6-7177d6d812cc
md"### snippet 2.6"

# ╔═╡ 3053171e-f094-11ea-171d-abc63f0f9766
p_grid = range(0, step=0.001, stop=1);      # Grid of 1001 steps

# ╔═╡ 75652c96-f08f-11ea-189b-d3e618f209a2
# all priors = 1.0

prior = ones(length(p_grid));

# ╔═╡ 7565dee8-f08f-11ea-00e4-ab58a09d11e3
# Binomial pdf

likelihood = [pdf(Binomial(9, p), 6) for p in p_grid];

# ╔═╡ 757a042c-f08f-11ea-32b7-8b10e91d9b9d
# Scale posterior such that they become probabilities
begin
	posterior = likelihood .* prior;
	posterior = posterior / sum(posterior);
end

# ╔═╡ 757ac0da-f08f-11ea-356e-51441c0702c6
N = 10000;         # Sample using the computed posterior values as weights

# ╔═╡ 758a9684-f08f-11ea-2d34-1f623f5f0774
samples = sample(p_grid, Weights(posterior), N);

# ╔═╡ 75946b00-f08f-11ea-3778-2d36aa151afb
chn = MCMCChains.Chains(reshape(samples, N, 1, 1), ["toss"]);

# ╔═╡ 7595cc16-f08f-11ea-21ef-1bfd716e8944
# Describe the chain

chn

# ╔═╡ 759ddcd0-f08f-11ea-079b-c11d8a66d53d
# Compute the MAP (maximum_a_posteriori) estimate

function loglik(x)
  ll = 0.0
  ll += log.(pdf.(Beta(1, 1), x[1]))
  ll += sum(log.(pdf.(Binomial(9, x[1]), repeat([6], 1))))
  -ll
end

# ╔═╡ 75a60f7e-f08f-11ea-1869-cf97d460cba1
opt = optimize(loglik, 0.0, 1.0)

# ╔═╡ 75b65508-f08f-11ea-15cd-7b01dde48e68
qmap = Optim.minimizer(opt)

# ╔═╡ 75b79198-f08f-11ea-056d-3f2c07bfe95d
# Show optimization results

opt

# ╔═╡ 75bf9262-f08f-11ea-3e65-ffde75e4b805
# Fit quadratic approcimation

quapfit = [qmap[1], std(samples, mean=qmap[1])]

# ╔═╡ 75c84510-f08f-11ea-27b6-57d559983738
p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 4);

# ╔═╡ 75ec0522-f08f-11ea-29c8-9fdca606cfd0
begin
	w = 6
	n = 9
	x = 0:0.01:1
end

# ╔═╡ 75ef0362-f08f-11ea-3645-e98e334241f0
p[1] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(-0.5, 1.0), 
  lab="Conjugate solution", leg=:topleft);

# ╔═╡ 75f9ea2a-f08f-11ea-232c-abb7e23fd457
density!(p[1], samples, lab="Sample density")

# ╔═╡ 7603dbca-f08f-11ea-3a8f-f3cd3f25ff35
# Distribution estimates copied from Turing quap()
# snippet 2.6

@model globethrowing(W, L) = begin
    p ~ Uniform(0, 1)
    W ~ Binomial(W + L, p)
end

# ╔═╡ 760e0958-f08f-11ea-12c1-bb2f7bc109b6
m = globethrowing(6, 3)

# ╔═╡ 76181db0-f08f-11ea-1afd-2de1941fb6d6
r = quap(m)

# ╔═╡ 762390aa-f08f-11ea-2fbb-cd5f4d92f5a4
d = Normal([r.coef.p][1], √collect(reshape(r.vcov, 1, 1))[1])

# ╔═╡ 762dee10-f08f-11ea-2fe2-5b2b4032a9b5
# quadratic approximation using Optim

p[2] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(-0.5, 1.0),
  lab="Conjugate solution", leg=:topleft);

# ╔═╡ 76387718-f08f-11ea-2eb6-23b4482d2623
plot!( p[2], x, pdf.(Normal( quapfit[1], quapfit[2] ) , x ),
  lab="Optim logpdf approx.");

# ╔═╡ 7642feae-f08f-11ea-13ae-e1d8b539fe80
plot!(p[2], x, pdf.(d, x), lab="Turing quap approx.")

# ╔═╡ 7658d45e-f08f-11ea-2f6b-b5f527ecace6
# quadratic approximation using StatisticalRethinking.jl quap()

df = DataFrame(:toss => samples);

# ╔═╡ 765da330-f08f-11ea-1543-0562390d470f
q = quap(df)

# ╔═╡ 76784712-f08f-11ea-012a-43cb94d9e255
p[3] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(-0.5, 1.0),
  lab="Conjugate solution", leg=:topleft);

# ╔═╡ 76846c72-f08f-11ea-1db9-1f60330ac3b2
plot!( p[3], x, pdf.(Normal(mean(q.toss), std(q.toss) ) , x ),
  lab="Stan quap approx.");

# ╔═╡ 7691174c-f08f-11ea-3472-fda8ccda44a5
plot!(p[3], x, pdf.(d, x), lab="Turing quap approx.")

# ╔═╡ 76aabf14-f08f-11ea-2f7a-c70588a9ac22
p[4] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(-0.5, 1.0),
  lab="Conjugate solution", leg=:topleft);

# ╔═╡ 76b7b098-f08f-11ea-031f-2d92bf266593
f = fit(Normal, samples)

# ╔═╡ 76c58112-f08f-11ea-1c3f-d542f6afaab9
plot!(p[4], x, pdf.(Normal( f.μ , f.σ ) , x ), lab="Normal MLE approx.")

# ╔═╡ 00cd4da6-f115-11ea-173b-2325e1576572
md"##### Summary plot Fig 2.7"

# ╔═╡ 76db71a0-f08f-11ea-1594-47a6175438eb
plot(p..., layout=(2, 2))

# ╔═╡ 76f96356-f08f-11ea-2b20-c3eb12968209
md"## End of clip-02-06-07t.jl"

# ╔═╡ Cell order:
# ╟─ca40913e-f093-11ea-28be-15a831f7c1c2
# ╠═a2e13c44-f08e-11ea-13df-53d884942d7e
# ╟─0c39592c-f094-11ea-3cb9-17977c11b058
# ╠═021cabc4-f08f-11ea-1fc2-f19729723f2c
# ╟─7554b2f8-f08f-11ea-38c6-7177d6d812cc
# ╠═3053171e-f094-11ea-171d-abc63f0f9766
# ╠═75652c96-f08f-11ea-189b-d3e618f209a2
# ╠═7565dee8-f08f-11ea-00e4-ab58a09d11e3
# ╠═757a042c-f08f-11ea-32b7-8b10e91d9b9d
# ╠═757ac0da-f08f-11ea-356e-51441c0702c6
# ╠═758a9684-f08f-11ea-2d34-1f623f5f0774
# ╠═75946b00-f08f-11ea-3778-2d36aa151afb
# ╠═7595cc16-f08f-11ea-21ef-1bfd716e8944
# ╠═759ddcd0-f08f-11ea-079b-c11d8a66d53d
# ╠═75a60f7e-f08f-11ea-1869-cf97d460cba1
# ╠═75b65508-f08f-11ea-15cd-7b01dde48e68
# ╠═75b79198-f08f-11ea-056d-3f2c07bfe95d
# ╠═75bf9262-f08f-11ea-3e65-ffde75e4b805
# ╠═75c84510-f08f-11ea-27b6-57d559983738
# ╠═75ec0522-f08f-11ea-29c8-9fdca606cfd0
# ╠═75ef0362-f08f-11ea-3645-e98e334241f0
# ╠═75f9ea2a-f08f-11ea-232c-abb7e23fd457
# ╠═7603dbca-f08f-11ea-3a8f-f3cd3f25ff35
# ╠═760e0958-f08f-11ea-12c1-bb2f7bc109b6
# ╠═76181db0-f08f-11ea-1afd-2de1941fb6d6
# ╠═762390aa-f08f-11ea-2fbb-cd5f4d92f5a4
# ╠═762dee10-f08f-11ea-2fe2-5b2b4032a9b5
# ╠═76387718-f08f-11ea-2eb6-23b4482d2623
# ╠═7642feae-f08f-11ea-13ae-e1d8b539fe80
# ╠═7658d45e-f08f-11ea-2f6b-b5f527ecace6
# ╠═765da330-f08f-11ea-1543-0562390d470f
# ╠═76784712-f08f-11ea-012a-43cb94d9e255
# ╠═76846c72-f08f-11ea-1db9-1f60330ac3b2
# ╠═7691174c-f08f-11ea-3472-fda8ccda44a5
# ╠═76aabf14-f08f-11ea-2f7a-c70588a9ac22
# ╠═76b7b098-f08f-11ea-031f-2d92bf266593
# ╠═76c58112-f08f-11ea-1c3f-d542f6afaab9
# ╟─00cd4da6-f115-11ea-173b-2325e1576572
# ╠═76db71a0-f08f-11ea-1594-47a6175438eb
# ╟─76f96356-f08f-11ea-2b20-c3eb12968209
