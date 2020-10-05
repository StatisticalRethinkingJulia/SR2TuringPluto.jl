### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ b80a0fc0-0675-11eb-1fcd-277d5ec1c289
using Pkg, DrWatson

# ╔═╡ b80a447c-0675-11eb-319c-01dc82d3fb7c
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ b221c3e6-0675-11eb-1ea1-4307128c6de3
md"## Clip-04-07-25t.jl"

# ╔═╡ b80ac79e-0675-11eb-2e70-17f8aa994dc0
md"## Clip-04-97-25t.jl"

# ╔═╡ b81d124e-0675-11eb-16e7-91a1a32f0369
md"### snippets 4.7 - 4.11"

# ╔═╡ b81df8d2-0675-11eb-33ad-819faf448b8a
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
end;

# ╔═╡ b82e9b88-0675-11eb-27b5-69985a392610
md"### snippets 4.12, 4.13"

# ╔═╡ b82f4f42-0675-11eb-3f77-d9e53da11c70
begin
	fig1 = plot(p -> pdf(Normal(178, 20), p), 100:250, leg=false, title="prior mu")
	fig2 = plot(p -> pdf(Uniform(0, 50), p), -10:60, leg=false, title="prior sigma")
	plot(fig1, fig2, layout=(1, 2))
end

# ╔═╡ b83efd16-0675-11eb-2e8d-fd33be380806
md"### snippet 4.14"

# ╔═╡ b84a549a-0675-11eb-08c9-b908f9639415
begin
	N = 1000
	d_μ = Normal(178, 20)
	sample_mu = rand(d_μ, N)
	d_σ = Uniform(0, 50)
	sample_sigma = rand(d_σ, N)
	prior_h = rand.(Normal.(sample_mu, sample_sigma))
	density(prior_h)
end

# ╔═╡ b85bf01a-0675-11eb-17e9-67c60998b96b
md"### snippet 4.15"

# ╔═╡ b8664768-0675-11eb-07de-555dbc1de52b
begin
	d_μ_100 = Normal(178, 100)
	sample_μ_100 = rand(d_μ_100, N)
end;

# ╔═╡ b86fd274-0675-11eb-2564-dfb336f2f5e6
begin
	sample_height_100 = rand.(Normal.(sample_μ_100, sample_sigma))
	density(sample_height_100)
end

# ╔═╡ 97fbe892-0676-11eb-0a85-6d793b380f4e
md"### snippets 4.16, 4.17"

# ╔═╡ 9b8b89fe-0676-11eb-1eaa-e986dcba4794
begin
	μ_grid_1 = range(150, 160, length = 50)
	σ_grid_1 = range(7, 9, length = 50)
	post_1 = [(; :μ => μ, :σ => σ) for μ in μ_grid_1, σ in σ_grid_1]  # this is `post`
	ll_1 = [loglikelihood(Normal(μ, σ), df.height) for (μ, σ) in post_1]
	prior_1 = [logpdf(d_μ, μ) + logpdf(d_σ, σ) for (μ, σ) in post_1]
	lposterior_1 = ll_1 .+ prior_1
	lposterior_1 .-= maximum(lposterior_1)
	posterior_1 = exp.(lposterior_1)
	wireframe(μ_grid_1, σ_grid_1, posterior_1)
end

# ╔═╡ b897b90a-0676-11eb-1fe9-cbd663edb4f0
md"### snippets 4.19 - 4.22"

# ╔═╡ c3a4dba2-0676-11eb-07fc-bfb72963914b
begin
	samples = sample(post_1, pweights(posterior_1), N) |> DataFrame
	scatter(samples.μ, samples.σ, alpha = 0.05, ms = 5, xlabel = "μ", ylabel = "σ")
end

# ╔═╡ fdacb630-0676-11eb-0217-f1e07c127aa2
begin
	fig3 = density(samples.μ, leg=false, title="samples.mu")
	fig4 = density(samples.σ, leg=false, title="samples.sigma")
	plot(fig3, fig4, layout=(1,2))
end

# ╔═╡ 643e622c-0677-11eb-187b-033fe474252d
quantile(samples.μ, (0.1, 0.9))

# ╔═╡ 6eed8892-0677-11eb-3c04-e9cb30ea2b66
quantile(samples.σ, (0.1, 0.9))

# ╔═╡ 776436b0-0677-11eb-19bf-bd953d980210
md"### snippets 4.23 - 4.25"

# ╔═╡ 9ba86ca0-0677-11eb-3304-118f06fda2f0
begin
	d3_height = df.height[1:20]
	μ_grid = range(140, 160, length = 50)
	σ_grid = range(4, 20, length = 50)
	post = [(μ = μ, σ = σ) for μ in μ_grid, σ in σ_grid]
	ll = [loglikelihood(Normal(μ, σ), d3_height) for (μ, σ) in post]
	prior = [logpdf(d_μ, μ) + logpdf(d_σ, σ) for (μ, σ) in post]
	lposterior = ll .+ prior
	lposterior .-= maximum(lposterior)
	posterior = exp.(lposterior)
	wireframe(μ_grid, σ_grid, posterior)
end

# ╔═╡ 8613457c-0675-11eb-10bd-217d514940ef
begin
	samples2 = sample(post, pweights(posterior), N) |> DataFrame
	scatter(samples2.μ, samples2.σ, alpha = 0.05, ms = 5)
end

# ╔═╡ d8ab7426-067a-11eb-0bac-1bb4747061b1
begin
	fig5 = density(samples2.μ, leg=false, title="samples2.mu")
	fig6 = density(samples2.σ, leg=false, title="samples2.sigma")
	plot(fig1, fig2, fig3, fig4, fig5, fig6, layout=(3,2))
end

# ╔═╡ d8d3d31c-067a-11eb-2a92-ef47acab4194
md"## End of clip-04-07-25t.jl"

# ╔═╡ Cell order:
# ╠═b221c3e6-0675-11eb-1ea1-4307128c6de3
# ╠═b80a0fc0-0675-11eb-1fcd-277d5ec1c289
# ╠═b80a447c-0675-11eb-319c-01dc82d3fb7c
# ╟─b80ac79e-0675-11eb-2e70-17f8aa994dc0
# ╟─b81d124e-0675-11eb-16e7-91a1a32f0369
# ╠═b81df8d2-0675-11eb-33ad-819faf448b8a
# ╟─b82e9b88-0675-11eb-27b5-69985a392610
# ╠═b82f4f42-0675-11eb-3f77-d9e53da11c70
# ╟─b83efd16-0675-11eb-2e8d-fd33be380806
# ╠═b84a549a-0675-11eb-08c9-b908f9639415
# ╟─b85bf01a-0675-11eb-17e9-67c60998b96b
# ╠═b8664768-0675-11eb-07de-555dbc1de52b
# ╠═b86fd274-0675-11eb-2564-dfb336f2f5e6
# ╟─97fbe892-0676-11eb-0a85-6d793b380f4e
# ╠═9b8b89fe-0676-11eb-1eaa-e986dcba4794
# ╟─b897b90a-0676-11eb-1fe9-cbd663edb4f0
# ╠═c3a4dba2-0676-11eb-07fc-bfb72963914b
# ╠═fdacb630-0676-11eb-0217-f1e07c127aa2
# ╠═643e622c-0677-11eb-187b-033fe474252d
# ╠═6eed8892-0677-11eb-3c04-e9cb30ea2b66
# ╟─776436b0-0677-11eb-19bf-bd953d980210
# ╠═9ba86ca0-0677-11eb-3304-118f06fda2f0
# ╠═8613457c-0675-11eb-10bd-217d514940ef
# ╠═d8ab7426-067a-11eb-0bac-1bb4747061b1
# ╟─d8d3d31c-067a-11eb-2a92-ef47acab4194
