### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ ea02e53e-e0cc-11ea-2265-3165c74e19b3
using Pkg, DrWatson

# ╔═╡ 02ea1c3e-e0cd-11ea-04d1-f34150f81c89
begin
	@quickactivate "StatisticalRethinkingTuring"
  	using Turing
	using StatisticalRethinking
end

# ╔═╡ d8ea5d90-e0cc-11ea-0d2d-25c807c1ae80
md"## Fig 2.5t"

# ╔═╡ 11533bb6-e0cd-11ea-331e-278de5d6b26f
md"### This clip is only intended to generate Fig 2.5."

# ╔═╡ 3067f2d0-e0cd-11ea-17d4-ab40276a0379
@model globe_toss(n, k) = begin
  theta ~ Beta(1, 1) # prior
  k ~ Binomial(n, theta) # model
  return k, theta
end;

# ╔═╡ a501be4a-fb4a-11ea-3820-ad09c47beb7a
begin
	p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 9)
	dens = Vector{DataFrame}(undef, 10)
	for n in 1:9
		p[n] = plot(xlims=(0.0, 1.0), ylims=(0.0, 3.0), leg=false)
		k = sum([1,0,1,1,1,0,1,0,1,1,0,1,1,1,0,1,0,1][1:n])
		chns = sample(globe_toss(n, k), NUTS(0.65), 1000)
		dfs=DataFrame(chns)
		if n == 1
			hline!([1.0], line=(:dash))
		else
			density!(dens[n][:, :theta], line=(:dash))
		end
		density!(dfs[:, :theta])
		dens[n+1] = dfs

	end
end


# ╔═╡ dd63b5f0-e0cd-11ea-0063-61ba73f99cac
plot(p..., layout=(3, 3))

# ╔═╡ ee6b3094-e0cd-11ea-1ceb-6f178f55cb23
md"## End of Fig2.5t.jl"

# ╔═╡ Cell order:
# ╟─d8ea5d90-e0cc-11ea-0d2d-25c807c1ae80
# ╠═ea02e53e-e0cc-11ea-2265-3165c74e19b3
# ╠═02ea1c3e-e0cd-11ea-04d1-f34150f81c89
# ╟─11533bb6-e0cd-11ea-331e-278de5d6b26f
# ╠═3067f2d0-e0cd-11ea-17d4-ab40276a0379
# ╠═a501be4a-fb4a-11ea-3820-ad09c47beb7a
# ╠═dd63b5f0-e0cd-11ea-0063-61ba73f99cac
# ╟─ee6b3094-e0cd-11ea-1ceb-6f178f55cb23
