### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 13280ba4-064a-11eb-04dd-5b2c0598d84a
using Pkg, DrWatson

# ╔═╡ 134f547a-064a-11eb-0750-b590037632f5
begin
	@quickactivate "StatisticalRethinkingTuring"
	using LinearAlgebra
	using Turing
	using StatisticalRethinking
	Turing.setprogress!(false)
end

# ╔═╡ b32d14e8-0648-11eb-15d4-7d13643a1531
md"## Clip-04-26-36t.jl"

# ╔═╡ 134fd792-064a-11eb-12dc-c5fe09197eba
md"### snippet 4.26"

# ╔═╡ 13549444-064a-11eb-08db-4fec5f434e9e
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
end;

# ╔═╡ 1363941c-064a-11eb-29c1-abbe84e92300
md"### snippet 4.27 - 4.29"

# ╔═╡ 13647ab2-064a-11eb-191c-73b11220323a
@model function ppl4_1(heights)
    μ ~ Normal(178, 20)
    σ ~ Uniform(0, 50)
    heights .~ Normal(μ, σ)
end

# ╔═╡ 137a6b92-064a-11eb-1b02-61006c9f3d8c
md"### snippet 4.30"

# ╔═╡ 136ffa54-064a-11eb-24a2-bbadb2cd9327
m4_1t = ppl4_1(df.height)

# ╔═╡ a452e82e-12fc-11eb-3f66-e96efc73e191
precis(m4_1t)

# ╔═╡ 1370b690-064a-11eb-01ff-ad478891bb93
md"##### Evaluating above cell might occasionally fail in Optim. Below cell should work."

# ╔═╡ 137b04a8-064a-11eb-3c7f-7959d865eda3
begin
	start = [mean(Normal(178, 20)), mean(Uniform(0, 50))]
	q4_1t = quap(m4_1t, start)
end

# ╔═╡ 13878f3e-064a-11eb-3d2d-833dab5fa561
md"### snippet 4.31"

# ╔═╡ 13886080-064a-11eb-014a-4badada265d5
@model function ppl4_2(heights)
    μ ~ Normal(178, 0.1)
    σ ~ Uniform(0, 50)
    heights .~ Normal(μ, σ)
end

# ╔═╡ 1394a660-064a-11eb-3b04-9fd3f1391d7f
begin
	m4_2t = ppl4_2(df.height)
	q4_2t = quap(m4_2t, NelderMead())
end

# ╔═╡ 3933a352-12fd-11eb-20bd-5bc9941a3fdc
typeof(q4_2t)

# ╔═╡ 5c64a180-12fd-11eb-2495-b7e35230241f
q4_1t.params

# ╔═╡ 1396d066-064a-11eb-16f6-41f3d6ed887e
md"### snippets 4.32, 4.33"

# ╔═╡ 13a1b544-064a-11eb-23eb-bd4cfd891d0f
q4_1t.vcov

# ╔═╡ 13a9efde-064a-11eb-00fe-99d6a9dae8c6
diag(q4_1t.vcov)

# ╔═╡ 78bb877e-064b-11eb-19e6-0f9da5cab17c
cov2cor(Matrix(q4_1t.vcov), sqrt.(diag(q4_1t.vcov)))

# ╔═╡ 13acf544-064a-11eb-2d81-3f22922e46ac
md"### snippets 4.34 - 4.36"

# ╔═╡ 13ba0c16-064a-11eb-0c6c-f19686f83ffd
begin
	res = rand(q4_1t.distr, 10_000)
	quap4_1t = DataFrame(res', ["μ", "σ"])
	Text(precis(quap4_1t; io=String))
end

# ╔═╡ 7cd372a2-12ff-11eb-1a05-bb7f9a95c1b7
begin
	chns4_1t = sample(m4_1t, NUTS(), 1000)
end

# ╔═╡ 13bbb6a6-064a-11eb-0658-99ca806b8a22
md"## End of clip-04-26-36t.jl"

# ╔═╡ Cell order:
# ╟─b32d14e8-0648-11eb-15d4-7d13643a1531
# ╠═13280ba4-064a-11eb-04dd-5b2c0598d84a
# ╠═134f547a-064a-11eb-0750-b590037632f5
# ╟─134fd792-064a-11eb-12dc-c5fe09197eba
# ╠═13549444-064a-11eb-08db-4fec5f434e9e
# ╟─1363941c-064a-11eb-29c1-abbe84e92300
# ╠═13647ab2-064a-11eb-191c-73b11220323a
# ╟─137a6b92-064a-11eb-1b02-61006c9f3d8c
# ╠═136ffa54-064a-11eb-24a2-bbadb2cd9327
# ╠═a452e82e-12fc-11eb-3f66-e96efc73e191
# ╟─1370b690-064a-11eb-01ff-ad478891bb93
# ╠═137b04a8-064a-11eb-3c7f-7959d865eda3
# ╟─13878f3e-064a-11eb-3d2d-833dab5fa561
# ╠═13886080-064a-11eb-014a-4badada265d5
# ╠═1394a660-064a-11eb-3b04-9fd3f1391d7f
# ╠═3933a352-12fd-11eb-20bd-5bc9941a3fdc
# ╠═5c64a180-12fd-11eb-2495-b7e35230241f
# ╟─1396d066-064a-11eb-16f6-41f3d6ed887e
# ╠═13a1b544-064a-11eb-23eb-bd4cfd891d0f
# ╠═13a9efde-064a-11eb-00fe-99d6a9dae8c6
# ╠═78bb877e-064b-11eb-19e6-0f9da5cab17c
# ╟─13acf544-064a-11eb-2d81-3f22922e46ac
# ╠═13ba0c16-064a-11eb-0c6c-f19686f83ffd
# ╠═7cd372a2-12ff-11eb-1a05-bb7f9a95c1b7
# ╟─13bbb6a6-064a-11eb-0658-99ca806b8a22
