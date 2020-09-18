### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 3c57ceb0-f934-11ea-1bef-c16b272d0c6b
using Pkg, DrWatson

# ╔═╡ 3c580a7e-f934-11ea-3b0c-4d5d7d5a7382
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ 83ad95c0-f933-11ea-1552-23ac7681803c
md"## Clip-02-03-06t.jl"

# ╔═╡ 3c588aa8-f934-11ea-0b0f-e32924dc7b0c
md"## snippet 2.3"

# ╔═╡ 3c65efac-f934-11ea-0719-b3745f3cd459
begin
	p_grid = range(0, 1, length = 20)
	prior = ones(20)
	likelihood = pdf.(Binomial.(9, p_grid), 6)
	posterior = likelihood .* prior
	posterior ./= sum(posterior)
end

# ╔═╡ 3c7fade8-f934-11ea-198a-85360398831c
md"## snippet 2.4"

# ╔═╡ 3c80d5f8-f934-11ea-03ed-6f4717473ad9
plot(p_grid, posterior, m = 3, legend = false,
    xlabel = "probability of water", ylabel = "posterior probability", title = "20 points")

# ╔═╡ 3c8eedf0-f934-11ea-28f8-47dca0c24ac7
md"## snippet 2.5"

# ╔═╡ 3c9939d6-f934-11ea-0ff1-5f6abcc5560e
begin
	prior2 = ifelse.(p_grid .< 0.5, 0, 1)
	prior3 = @. exp(-5 * abs(p_grid - 0.5))  # @. means broadcast everything that follows
end

# ╔═╡ 3c9b1db4-f934-11ea-3d97-f51c2acc9037
md"## snippet 2.6"

# ╔═╡ 3ca3ac90-f934-11ea-1a8e-edfa2c664f1a
@model globethrowing(W, L) = begin
    p ~ Uniform(0, 1)
    W ~ Binomial(W + L, p)
end

# ╔═╡ 3caf45dc-f934-11ea-2d1d-796707d89526
m = globethrowing(6, 3)

# ╔═╡ 3cb18732-f934-11ea-3ebf-991c7b602d12
r = quap(m)

# ╔═╡ 3cb8986c-f934-11ea-1ffe-693fa1640736
d = MvNormal([r.coef.p], √collect(reshape(r.vcov, 1, 1)))

# ╔═╡ 3cbf6b30-f934-11ea-3b10-433589e36648
s = rand(d, 10000)';

# ╔═╡ 3cc6c8d8-f934-11ea-0a23-8db011bce2bf
histogram(collect(s), normalize = :probability)

# ╔═╡ 3cd586b6-f934-11ea-3820-b9ff7509f935
plot(p_grid, posterior, m = 3)

# ╔═╡ 3ce617ae-f934-11ea-0e4a-ff1ec9ee2610
md"## End clip-02-03-06t.jl"

# ╔═╡ Cell order:
# ╟─83ad95c0-f933-11ea-1552-23ac7681803c
# ╠═3c57ceb0-f934-11ea-1bef-c16b272d0c6b
# ╠═3c580a7e-f934-11ea-3b0c-4d5d7d5a7382
# ╟─3c588aa8-f934-11ea-0b0f-e32924dc7b0c
# ╠═3c65efac-f934-11ea-0719-b3745f3cd459
# ╟─3c7fade8-f934-11ea-198a-85360398831c
# ╠═3c80d5f8-f934-11ea-03ed-6f4717473ad9
# ╟─3c8eedf0-f934-11ea-28f8-47dca0c24ac7
# ╠═3c9939d6-f934-11ea-0ff1-5f6abcc5560e
# ╟─3c9b1db4-f934-11ea-3d97-f51c2acc9037
# ╠═3ca3ac90-f934-11ea-1a8e-edfa2c664f1a
# ╠═3caf45dc-f934-11ea-2d1d-796707d89526
# ╠═3cb18732-f934-11ea-3ebf-991c7b602d12
# ╠═3cb8986c-f934-11ea-1ffe-693fa1640736
# ╠═3cbf6b30-f934-11ea-3b10-433589e36648
# ╠═3cc6c8d8-f934-11ea-0a23-8db011bce2bf
# ╠═3cd586b6-f934-11ea-3820-b9ff7509f935
# ╟─3ce617ae-f934-11ea-0e4a-ff1ec9ee2610
