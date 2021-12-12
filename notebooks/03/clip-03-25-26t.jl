### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 4bba7820-f52f-11ea-2ffd-879f7444eb9c
using Pkg, DrWatson

# ╔═╡ 4bbabd56-f52f-11ea-0da9-55fcb93b0375
begin
	@quickactivate "SR2TuringPluto"
	using StatisticalRethinking
end

# ╔═╡ 0622be1c-f52f-11ea-026d-39e1d2e6b03c
md"## Clip-03-25-26t.jl"

# ╔═╡ 4bbb41ba-f52f-11ea-18bd-039df6e7d6cb
md"### snippet 3.25"

# ╔═╡ 4bc8abfc-f52f-11ea-1a7e-515b267bb1e1
begin
	N = 10000
	w = rand(Binomial(9, 0.6), N);
	h1 = histogram(w; normalize=:probability, leg=false)
end

# ╔═╡ 4bc93f7c-f52f-11ea-10c7-1595439f2b94
md"##### Generate the samples."

# ╔═╡ 4bd64e92-f52f-11ea-3326-fb4b8eb664a4
begin
	p_grid = range(0, step=0.001, stop=1)
	prior = ones(length(p_grid))
	likelihood = [pdf(Binomial(9, p), 6) for p in p_grid]
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
	samples = sample(p_grid, Weights(posterior), N)
end;

# ╔═╡ 4bd71c82-f52f-11ea-1706-876cbc5db9d5
md"### snippet 3.26"

# ╔═╡ 4bdfebf0-f52f-11ea-22df-43be326fcd25
begin
	d = rand.(Binomial.(9, samples));
	h2 = histogram(d; normalize=:probability, 
	  bins=-0.5:1:9.5, leg=false, xticks=0:9, bar_width=0.2)
end

# ╔═╡ 4be073c2-f52f-11ea-2ad1-5b555bbb69be
plot(h1, h2, layout=(1, 2))

# ╔═╡ 4bec4ae4-f52f-11ea-34c1-cdbc9b668632
md"## End of clip-03-25-26t.jl"

# ╔═╡ Cell order:
# ╟─0622be1c-f52f-11ea-026d-39e1d2e6b03c
# ╠═4bba7820-f52f-11ea-2ffd-879f7444eb9c
# ╠═4bbabd56-f52f-11ea-0da9-55fcb93b0375
# ╟─4bbb41ba-f52f-11ea-18bd-039df6e7d6cb
# ╠═4bc8abfc-f52f-11ea-1a7e-515b267bb1e1
# ╟─4bc93f7c-f52f-11ea-10c7-1595439f2b94
# ╠═4bd64e92-f52f-11ea-3326-fb4b8eb664a4
# ╟─4bd71c82-f52f-11ea-1706-876cbc5db9d5
# ╠═4bdfebf0-f52f-11ea-22df-43be326fcd25
# ╠═4be073c2-f52f-11ea-2ad1-5b555bbb69be
# ╟─4bec4ae4-f52f-11ea-34c1-cdbc9b668632
