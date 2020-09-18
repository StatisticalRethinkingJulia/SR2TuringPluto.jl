### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 9e7c03f0-f771-11ea-37c6-7bb5bdfee6b3
using Pkg, DrWatson

# ╔═╡ 9e7c4612-f771-11ea-18ab-53b0327f84f9
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ da886344-f770-11ea-2dce-c32c5b634485
md"## Fig2.9t.jl"

# ╔═╡ 730ca458-f772-11ea-37b7-3b1c57b1cf3c
md"##### Comparisons between different quadratic approximations."

# ╔═╡ 9e7cd226-f771-11ea-16ff-5f359383d3df
md"##### Grid of 1001 steps."

# ╔═╡ 9e8b4d24-f771-11ea-13d9-53d4d897b6f9
p_grid = range(0, step=0.001, stop=1);

# ╔═╡ 9e8bebee-f771-11ea-19bb-45a7a7d0d1c7
md"##### Set all priors = 1.0."

# ╔═╡ 9e9831d8-f771-11ea-1f51-0f8f67412b1a
prior = ones(length(p_grid));

# ╔═╡ 9e98c88c-f771-11ea-1b71-6954d845ece0
md"##### Binomial pdf."

# ╔═╡ 9ea43208-f771-11ea-0928-596e88c658ae
likelihood = pdf.(Binomial.(9, p_grid), 6);

# ╔═╡ 9ea4dfbe-f771-11ea-38dc-9f3b1e72cf4c
md"##### As Uniform prior has been used, unstandardized posterior is equal to likelihood."

# ╔═╡ 9eb26684-f771-11ea-2f72-9158a703be12
begin
	posterior = likelihood .* prior
	
	# Scale posterior such that they become probabilities."
	
	posterior = posterior / sum(posterior)
end;

# ╔═╡ 9eca0c1a-f771-11ea-267d-519ccb35777e
md"##### Sample using the computed posterior values as weights."

# ╔═╡ 9ed22960-f771-11ea-03aa-dfe49b859ac9
begin
	N = 10000
	samples = sample(p_grid, Weights(posterior), N);
	chn = MCMCChains.Chains(reshape(samples, N, 1, 1), ["toss"]);
end;

# ╔═╡ 9eda9df2-f771-11ea-2b33-0d04da6b3b8e
md"##### Describe the chain."

# ╔═╡ 9ee2b730-f771-11ea-0cb3-4380b1093306
chn

# ╔═╡ 9ee425c0-f771-11ea-1eec-2966da179664
md"##### Compute the MAP (maximum a posteriori) estimate."

# ╔═╡ 9eebf08e-f771-11ea-0b37-e39eae2dfac4
function loglik(x)
  ll = 0.0
  ll += log.(pdf.(Beta(1, 1), x[1]))
  ll += sum(log.(pdf.(Binomial(9, x[1]), repeat([6], 1))))
  -ll
end

# ╔═╡ 9ef4614c-f771-11ea-1e8f-1d0aa4f1abdf
opt = optimize(loglik, 0.0, 1.0)

# ╔═╡ 9efcc1ac-f771-11ea-0d80-f96d420ba180
qmap = Optim.minimizer(opt)

# ╔═╡ 9f1c1638-f771-11ea-3c94-b7c61d2d946d
md"##### Fit Optim quadratic approximation."

# ╔═╡ 9f1fd692-f771-11ea-20d4-1f85875e38e3
quapfit = [qmap[1], std(samples, mean=qmap[1])]

# ╔═╡ 9f2aa07c-f771-11ea-1725-1134b6105e6a
f = Vector{Plots.Plot{Plots.GRBackend}}(undef, 4);

# ╔═╡ bd367b42-f9ef-11ea-17a8-3f634500102e
md"##### Fit Turing quap approximation."

# ╔═╡ adf1a79e-f9ed-11ea-385a-e3a66c3c0901
begin
	w = 6
	l = 3
	n = 9
	x = 0:0.01:1
end

# ╔═╡ b4c7fc36-f9ec-11ea-1fff-49f9e92ceab4
@model globethrowing(w, l) = begin
    p ~ Uniform(0, 1)
    w ~ Binomial(w + l, p)
end

# ╔═╡ b41c2c12-f9ec-11ea-2c5c-f1e302d70c3a
begin
	m = globethrowing(6, 3)
	r = quap(m)
end

# ╔═╡ 9f3619ca-f771-11ea-24c4-39ec211debfc
md"##### Setup the plots."

# ╔═╡ 9f405a8e-f771-11ea-2df4-914745947f38
begin
	f[1] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(-0.5, 1.0), 
	  lab="Conjugate solution", leg=:topleft)
	density!(f[1], samples, lab="Sample density")

	# quadratic approximation using Optim

	f[2] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(-0.5, 1.0),
	  lab="Conjugate solution", leg=:topleft)
	plot!( f[2], x, pdf.(Normal( quapfit[1], quapfit[2] ) , x ),
	  lab="Optim logpdf approx.")

	# Quadratic approximation using Turing.jl quap()

	f[3] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(-0.5, 1.0),
	  lab="Conjugate solution", leg=:topleft)
	plot!( f[3], x, pdf.(Normal([r.coef.p][1], √collect(reshape(r.vcov, 1, 1))[1]) , x ), 
		lab="Turing quap approx.")

	# MLE approximation

	f[4] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(-0.5, 1.0),
	  lab="Conjugate solution", leg=:topleft)
	fi = fit(Normal, samples)
	plot!(f[4], x, pdf.(Normal( fi.μ , fi.σ ) , x ), lab="Normal MLE approx.")
end;

# ╔═╡ 9f4b6bae-f771-11ea-0e74-832e38739a30
plot(f..., layout=(2, 2))

# ╔═╡ 9f55638e-f771-11ea-348b-1b39f9b07e41
md"## End of Fig2.9t.jl"

# ╔═╡ Cell order:
# ╟─da886344-f770-11ea-2dce-c32c5b634485
# ╟─730ca458-f772-11ea-37b7-3b1c57b1cf3c
# ╠═9e7c03f0-f771-11ea-37c6-7bb5bdfee6b3
# ╠═9e7c4612-f771-11ea-18ab-53b0327f84f9
# ╠═9e7cd226-f771-11ea-16ff-5f359383d3df
# ╠═9e8b4d24-f771-11ea-13d9-53d4d897b6f9
# ╟─9e8bebee-f771-11ea-19bb-45a7a7d0d1c7
# ╠═9e9831d8-f771-11ea-1f51-0f8f67412b1a
# ╟─9e98c88c-f771-11ea-1b71-6954d845ece0
# ╠═9ea43208-f771-11ea-0928-596e88c658ae
# ╟─9ea4dfbe-f771-11ea-38dc-9f3b1e72cf4c
# ╠═9eb26684-f771-11ea-2f72-9158a703be12
# ╟─9eca0c1a-f771-11ea-267d-519ccb35777e
# ╠═9ed22960-f771-11ea-03aa-dfe49b859ac9
# ╟─9eda9df2-f771-11ea-2b33-0d04da6b3b8e
# ╠═9ee2b730-f771-11ea-0cb3-4380b1093306
# ╠═9ee425c0-f771-11ea-1eec-2966da179664
# ╠═9eebf08e-f771-11ea-0b37-e39eae2dfac4
# ╠═9ef4614c-f771-11ea-1e8f-1d0aa4f1abdf
# ╠═9efcc1ac-f771-11ea-0d80-f96d420ba180
# ╟─9f1c1638-f771-11ea-3c94-b7c61d2d946d
# ╠═9f1fd692-f771-11ea-20d4-1f85875e38e3
# ╠═9f2aa07c-f771-11ea-1725-1134b6105e6a
# ╟─bd367b42-f9ef-11ea-17a8-3f634500102e
# ╠═adf1a79e-f9ed-11ea-385a-e3a66c3c0901
# ╠═b4c7fc36-f9ec-11ea-1fff-49f9e92ceab4
# ╠═b41c2c12-f9ec-11ea-2c5c-f1e302d70c3a
# ╟─9f3619ca-f771-11ea-24c4-39ec211debfc
# ╠═9f405a8e-f771-11ea-2df4-914745947f38
# ╠═9f4b6bae-f771-11ea-0e74-832e38739a30
# ╟─9f55638e-f771-11ea-348b-1b39f9b07e41
