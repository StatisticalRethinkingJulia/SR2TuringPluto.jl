### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 9e6cead1-18cc-4ff1-bb4a-66e3ca8f37da
using Pkg

# ╔═╡ df0344f8-dcec-42fe-869d-2e89369d2201
#Pkg.activate(expanduser("~/.julia/dev/SR2TuringPluto"))

# ╔═╡ 39bd8de1-6664-4ffe-abcc-10b9f4788f36
begin
	using Distributions
	using Optim
	using Turing
	using Logging
	using LaTeXStrings
	using StatisticalRethinking
	using StatisticalRethinkingPlots
end

# ╔═╡ 6a78dd3a-9770-40d0-b344-e7678801499c
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

# ╔═╡ c1cec086-0ca1-449a-b9d6-e1366d09f811
begin
	Logging.disable_logging(Logging.Warn)
end;

# ╔═╡ 5826ebdc-3ecf-48ff-ae46-91ac572db0b3
md"### Code 2.1"

# ╔═╡ 26050f29-8713-48bb-a279-abd42cabb959
begin
	ways = [0, 3, 8, 9, 0]
	ways = ways ./ sum(ways)
end

# ╔═╡ 54e330b9-814a-4b34-9486-9d7f58019219
md"### Code 2.2"

# ╔═╡ 663bb843-fc37-4242-8df6-3099bc0d0d56
let
	b = Binomial(9, 0.5)
	pdf(b, 6)
end

# ╔═╡ 3e2fbbbc-39ce-405c-ad5c-6e3cc7c23611
md"### Codes 2.3 & 2.4"

# ╔═╡ f8994f78-08d5-4745-ad85-686a33aa6c57
md"#### Grid, prior, likelihood and posterior."

# ╔═╡ f8146fd4-586a-4b70-8d68-f96bb3a889c9
let
	size = 20
	p_grid = range(0, 1; length=size)
	prior = ones(size)
	likelyhood = [pdf(Binomial(9, p), 6) for p in p_grid]
	unstd_posterior = likelyhood .* prior
	posterior = unstd_posterior / sum(unstd_posterior)

	plot(p_grid, posterior; 
		xlabel="probability of water", 
		ylabel="posterior probability",
		title="$size points",
		markershape=:circle)
end

# ╔═╡ 6c52cc82-dae3-48d3-ab55-dd0991285600
md"#### Code 2.5"

# ╔═╡ 61dac78d-e092-413b-9044-3078d35a1abf
let
	size = 20
	p_grid = range(0, 1; length=size)
	prior = convert(Vector{AbstractFloat}, p_grid .>= 0.5)
	likelyhood = [pdf(Binomial(9, p), 6) for p in p_grid]
	unstd_posterior = likelyhood .* prior
	posterior = unstd_posterior / sum(unstd_posterior)

	plot(p_grid, posterior; 
		xlabel="probability of water", 
		ylabel="posterior probability",
		title="$size points",
		markershape=:circle)
end

# ╔═╡ 212cd629-9e1e-41b0-b28e-71f1c99e93f6
md"#### Another prior to try."

# ╔═╡ 2d6bacd4-6a5e-4857-accc-23dda28b69e5
let
	size = 20
	p_grid = range(0, 1; length=size)
	prior = exp.(-5*abs.(p_grid .- 0.5))
	likelyhood = [pdf(Binomial(9, p), 6) for p in p_grid]
	unstd_posterior = likelyhood .* prior
	posterior = unstd_posterior / sum(unstd_posterior)

	plot(p_grid, posterior; 
		xlabel="probability of water", 
		ylabel="posterior probability",
		title="$size points",
		markershape=:circle)
end

# ╔═╡ 14940af1-2cbb-4527-8a0c-83a54e305aac
md"### Code 2.6"

# ╔═╡ c34c3d5a-f28d-4819-92d7-53103fe22be7
md"### Code 2.7"

# ╔═╡ c39d8158-d6fe-4cc3-9ba2-7a418815341a
md"### Analytical calculation."

# ╔═╡ 443044d5-02a3-4bfc-955b-54b0c12c63e2
begin
	W = 6
	L = 3
	x = collect(range(0, 1; length=101))
	b1 = Beta(W+1, L+1)
	plot(x, pdf.(b1, x); label = "β")
end

# ╔═╡ 3276b558-7dc7-4869-a083-472716a8a2ab
md"### Quadratic approximation."

# ╔═╡ 20e03a5e-3321-4019-97a5-44a0fd77a831
md"### Code 2.8"

# ╔═╡ 4a47f055-842c-451c-bc75-072708588f9a
begin
	n_samples = 1000
	p = Vector{Float64}(undef, n_samples)
	p[1] = 0.5

	for i ∈ 2:n_samples
		p_old = p[i-1]
		p_new = rand(Normal(p_old, 0.1))
		if p_new < 0
			p_new = abs(p_new)
		elseif p_new > 1
			p_new = 2-p_new
		end

		q0 = pdf(Binomial(W+L, p_old), W)
		q1 = pdf(Binomial(W+L, p_new), W)
		u = rand(Uniform())
		p[i] = (u < q1 / q0) ? p_new : p_old
	end
end

# ╔═╡ 7b0a5745-d8f8-4475-bd68-8ee662a61fc1
@model function m2_0(W, L)
    p ~ Uniform(0, 1)
    W ~ Binomial(W + L, p)
end

# ╔═╡ e15efdfa-e64f-4f5b-a94c-6c5e135a9ec2
let
	map2_0 = optimize(m2_0(W, L), MAP())
	μ = coef(map2_0)[:p]
	σ = sqrt(vcov(map2_0)[:p, :p])
	b2 = Normal(μ, σ)
	plot!(x, pdf.(b2, x); style=:dash,
		label="Normal($(round(μ; digits=2)), $(round(σ; digits=2)))")
end

# ╔═╡ 8c7de4e1-b5b1-44c4-9718-b8044a6397ef
md"### Code 2.9"

# ╔═╡ 4e7745ef-50f5-4437-9223-cec855a3d2f1
begin
	density(p; label = "MCMC")
	b3 = Beta(W+1, L+1)
	plot!(x, pdf.(b3, x); label = L"\mathcal{Beta}", style=:dash)
end

# ╔═╡ a06463e3-4d9a-45cd-a0bc-3e44d37ad73c
let
	n = 9
	p = 0.66
	μ = n * p
	σ = sqrt(n * p * ( 1 - p ))
	(μ = μ, σ = σ)
end

# ╔═╡ Cell order:
# ╠═6a78dd3a-9770-40d0-b344-e7678801499c
# ╠═9e6cead1-18cc-4ff1-bb4a-66e3ca8f37da
# ╠═df0344f8-dcec-42fe-869d-2e89369d2201
# ╠═39bd8de1-6664-4ffe-abcc-10b9f4788f36
# ╠═c1cec086-0ca1-449a-b9d6-e1366d09f811
# ╟─5826ebdc-3ecf-48ff-ae46-91ac572db0b3
# ╠═26050f29-8713-48bb-a279-abd42cabb959
# ╟─54e330b9-814a-4b34-9486-9d7f58019219
# ╠═663bb843-fc37-4242-8df6-3099bc0d0d56
# ╟─3e2fbbbc-39ce-405c-ad5c-6e3cc7c23611
# ╟─f8994f78-08d5-4745-ad85-686a33aa6c57
# ╠═f8146fd4-586a-4b70-8d68-f96bb3a889c9
# ╟─6c52cc82-dae3-48d3-ab55-dd0991285600
# ╠═61dac78d-e092-413b-9044-3078d35a1abf
# ╠═212cd629-9e1e-41b0-b28e-71f1c99e93f6
# ╠═2d6bacd4-6a5e-4857-accc-23dda28b69e5
# ╟─14940af1-2cbb-4527-8a0c-83a54e305aac
# ╠═7b0a5745-d8f8-4475-bd68-8ee662a61fc1
# ╟─c34c3d5a-f28d-4819-92d7-53103fe22be7
# ╟─c39d8158-d6fe-4cc3-9ba2-7a418815341a
# ╠═443044d5-02a3-4bfc-955b-54b0c12c63e2
# ╟─3276b558-7dc7-4869-a083-472716a8a2ab
# ╠═e15efdfa-e64f-4f5b-a94c-6c5e135a9ec2
# ╟─20e03a5e-3321-4019-97a5-44a0fd77a831
# ╠═4a47f055-842c-451c-bc75-072708588f9a
# ╟─8c7de4e1-b5b1-44c4-9718-b8044a6397ef
# ╠═4e7745ef-50f5-4437-9223-cec855a3d2f1
# ╠═a06463e3-4d9a-45cd-a0bc-3e44d37ad73c
