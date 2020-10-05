
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Clip-04-07-25t.jl"

md"## Clip-04-97-25t.jl"

md"### snippets 4.7 - 4.11"

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
end;

md"### snippets 4.12, 4.13"

begin
	fig1 = plot(p -> pdf(Normal(178, 20), p), 100:250, leg=false, title="prior mu")
	fig2 = plot(p -> pdf(Uniform(0, 50), p), -10:60, leg=false, title="prior sigma")
	plot(fig1, fig2, layout=(1, 2))
end

md"### snippet 4.14"

begin
	N = 1000
	d_μ = Normal(178, 20)
	sample_mu = rand(d_μ, N)
	d_σ = Uniform(0, 50)
	sample_sigma = rand(d_σ, N)
	prior_h = rand.(Normal.(sample_mu, sample_sigma))
	density(prior_h)
end

md"### snippet 4.15"

begin
	d_μ_100 = Normal(178, 100)
	sample_μ_100 = rand(d_μ_100, N)
end;

begin
	sample_height_100 = rand.(Normal.(sample_μ_100, sample_sigma))
	density(sample_height_100)
end

md"### snippets 4.16, 4.17"

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

md"### snippets 4.19 - 4.22"

begin
	samples = sample(post_1, pweights(posterior_1), N) |> DataFrame
	scatter(samples.μ, samples.σ, alpha = 0.05, ms = 5, xlabel = "μ", ylabel = "σ")
end

begin
	fig3 = density(samples.μ, leg=false, title="samples.mu")
	fig4 = density(samples.σ, leg=false, title="samples.sigma")
	plot(fig3, fig4, layout=(1,2))
end

quantile(samples.μ, (0.1, 0.9))

quantile(samples.σ, (0.1, 0.9))

md"### snippets 4.23 - 4.25"

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

begin
	samples2 = sample(post, pweights(posterior), N) |> DataFrame
	scatter(samples2.μ, samples2.σ, alpha = 0.05, ms = 5)
end

begin
	fig5 = density(samples2.μ, leg=false, title="samples2.mu")
	fig6 = density(samples2.σ, leg=false, title="samples2.sigma")
	plot(fig1, fig2, fig3, fig4, fig5, fig6, layout=(3,2))
end

md"## End of clip-04-07-25t.jl"

