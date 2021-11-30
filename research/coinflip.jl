### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# ╔═╡ 73094af4-3595-11eb-269d-e3bf33913758
using Pkg, DrWatson

# ╔═╡ 73323464-3595-11eb-3be0-b9a9f1ea2f8f
begin
    @quickactivate "SR2TuringPluto"
    using Turing
    using StatisticalRethinking
end

# ╔═╡ 5550278a-3595-11eb-0e3e-f12ab620d448
md"## Coinflip.jl"

# ╔═╡ 7332c5fa-3595-11eb-288b-dd932e0953bf
struct MyDist{T} <: ContinuousUnivariateDistribution
    θ::T
end

# ╔═╡ 734cff6a-359a-11eb-2cde-4756a76d07a7
md"##### I suspect that the for loop in coinflip1 performs a separate evaluation on each iteration with overhead. Perhaps there is someway to improve performance on the Turing side. As a workaround, you can create a new distribution type and place the for loop inside the logpdf function. I use this in cases where
1. I am defining a custom distribution,
2. I am using a non-standard data type,
3. or a for loop is a more natural way to define the model."

# ╔═╡ cf9e35ae-359a-11eb-02ec-0550744aed96
md"##### Below some examples for comparison."

# ╔═╡ 73435406-3595-11eb-3e69-9946981f886d
begin
	import Distributions: logpdf	
	function logpdf(d::MyDist, data::Array{Bool,1})
		LL = 0.0
		for i in 1:length(data)
			LL += logpdf(Bernoulli(d.θ), data[i])
		end
		return LL
	end
end

# ╔═╡ 73450f44-3595-11eb-2707-65c5f3ff8451
begin
	import Distributions: loglikelihood
	loglikelihood(d::MyDist, data::Array{Bool,1}) = logpdf(d, data)
end

# ╔═╡ 7350c640-3595-11eb-1f8c-41159f4490cf
begin
    p_true = 0.5
    n = 1000
    Random.seed!(12)
    data = rand(Bernoulli(p_true), n)
end;

# ╔═╡ 73515718-3595-11eb-0c71-b3c8fae4e95b
begin
    n_samples = 2000
    delta = .65
    n_adapt = 1000
    config = NUTS(n_adapt, delta)
end;

# ╔═╡ 735effce-3595-11eb-06c4-9b8202849af3
@model function coinflip1(y)
    p ~ Beta(1.,1.)
    for i in 1:length(y)
        y[i] ~ Bernoulli(p)
    end
end

# ╔═╡ 735fa414-3595-11eb-3c38-8d1ba76fa28e
md"##### Is coinflip2 code correct (no documentation, but it works)?"

# ╔═╡ 736e5b4c-3595-11eb-1ff2-5b2ca39213ac
@model function coinflip2(y)
    p ~ Beta(1.,1.)
    y ~ Bernoulli(p)
end

# ╔═╡ 737a61da-3595-11eb-3f65-ab61047f1fbc
md"##### This model is from 'DynamicPPL: Stan-like Speed for Dynamic Probabilistic Models'"

# ╔═╡ 737d45c6-3595-11eb-0fec-3b044ca5e64c
@model function coinflip3(y)
    p ~ Beta(1.,1.)
    y .~ Bernoulli(p)
end

# ╔═╡ 7b19ca74-3599-11eb-1370-1f773dc7ae9d
md"##### This code uses above defined definition."

# ╔═╡ 738d7f04-3595-11eb-1b38-0315d5ae6083
@model function coinflip4(y)
    p ~ Beta(1.,1.)
    y ~ MyDist(p)
end

# ╔═╡ 4db69812-3598-11eb-31ea-a152bc58aca9
md"#### Note: Output of `@time ...` to Terminal (or in  right lower cell corner)."

# ╔═╡ 73944a32-3595-11eb-21a4-4f76a53efac2
@time chain1 = sample(coinflip1(data), config, n_samples);

# ╔═╡ 49bdbc14-359c-11eb-3e77-fb604ac0bf4f
@time chain2 = sample(coinflip2(data), config, n_samples);

# ╔═╡ 49bde4d2-359c-11eb-3440-0b244f8f784d
@time chain3 = sample(coinflip3(data), config, n_samples);

# ╔═╡ 49be84fa-359c-11eb-0f24-37ff390dc0b4
@time chain4 = sample(coinflip4(data), config, n_samples);

# ╔═╡ 90e7082e-3598-11eb-2ecf-f78e174b3dc8
begin
	chns = chainscat(chain1, chain2, chain3, chain4)
	CHNS(chns)
end

# ╔═╡ 73bbc712-3595-11eb-1eab-5fb16c871e73
md"## End of coinflip.jl"

# ╔═╡ Cell order:
# ╟─5550278a-3595-11eb-0e3e-f12ab620d448
# ╠═73094af4-3595-11eb-269d-e3bf33913758
# ╠═73323464-3595-11eb-3be0-b9a9f1ea2f8f
# ╠═7332c5fa-3595-11eb-288b-dd932e0953bf
# ╟─734cff6a-359a-11eb-2cde-4756a76d07a7
# ╟─cf9e35ae-359a-11eb-02ec-0550744aed96
# ╠═73435406-3595-11eb-3e69-9946981f886d
# ╠═73450f44-3595-11eb-2707-65c5f3ff8451
# ╠═7350c640-3595-11eb-1f8c-41159f4490cf
# ╠═73515718-3595-11eb-0c71-b3c8fae4e95b
# ╠═735effce-3595-11eb-06c4-9b8202849af3
# ╟─735fa414-3595-11eb-3c38-8d1ba76fa28e
# ╠═736e5b4c-3595-11eb-1ff2-5b2ca39213ac
# ╟─737a61da-3595-11eb-3f65-ab61047f1fbc
# ╠═737d45c6-3595-11eb-0fec-3b044ca5e64c
# ╟─7b19ca74-3599-11eb-1370-1f773dc7ae9d
# ╠═738d7f04-3595-11eb-1b38-0315d5ae6083
# ╟─4db69812-3598-11eb-31ea-a152bc58aca9
# ╠═73944a32-3595-11eb-21a4-4f76a53efac2
# ╠═49bdbc14-359c-11eb-3e77-fb604ac0bf4f
# ╠═49bde4d2-359c-11eb-3440-0b244f8f784d
# ╠═49be84fa-359c-11eb-0f24-37ff390dc0b4
# ╠═90e7082e-3598-11eb-2ecf-f78e174b3dc8
# ╟─73bbc712-3595-11eb-1eab-5fb16c871e73
