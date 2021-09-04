### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ e5368ab2-eb97-11ea-1755-f5b9bae4e5c3
using Pkg, DrWatson

# ╔═╡ cb4a13ee-eb97-11ea-373e-474e38399b44
begin
	@quickactivate "StatisticalRethinkingTuring"
	using GLM
	using StatisticalRethinking
end

# ╔═╡ 8b40ccde-eb97-11ea-0553-f5def6476477
md"## Clip-00-04-05t.jl"

# ╔═╡ 3d8845e0-eb98-11ea-128a-93226e62b034
md"### snippet 0.5"

# ╔═╡ 9b9bf518-eb97-11ea-2876-bd1758c8afd1
md"##### Snippet 0.5 is replaced by below lines."

# ╔═╡ fe16f312-eb97-11ea-1272-fb0eb517fbe2
md"### snippet 0.4"

# ╔═╡ 1059f184-eb98-11ea-1108-a50008b9be0c
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame; delim=';')
	df = filter(row -> row[:age] >= 18, df);
end;

# ╔═╡ eb567d8c-122d-11eb-0f13-f7ea2dc1f805
Text(precis(df; io=String))

# ╔═╡ 37924058-eb98-11ea-3a78-390939122024
md"##### Fit a linear regression of distance on speed."

# ╔═╡ 61681d9c-eb98-11ea-1411-df9bdf7296ae
m = lm(@formula(height ~ weight), df)

# ╔═╡ 686ef54a-eb98-11ea-1b68-d33aa9873780
md"##### Estimated coefficients from the model."

# ╔═╡ 5b284ec8-eb97-11ea-142f-85ff5053bd2b
coef(m)

# ╔═╡ b0943d46-eb98-11ea-1a87-8d8b461219a8
md"##### Plot residuals against height."

# ╔═╡ d7308a84-eb98-11ea-2299-95175b5cf27d
scatter( df.height, residuals(m), xlab="Height",
  ylab="Model residual values", lab="Model residuals", leg=:bottomright)

# ╔═╡ e17084e2-eb98-11ea-34ef-b1ebfab71041
md"## End of clip-00-04-05t.jl"

# ╔═╡ Cell order:
# ╟─8b40ccde-eb97-11ea-0553-f5def6476477
# ╟─3d8845e0-eb98-11ea-128a-93226e62b034
# ╟─9b9bf518-eb97-11ea-2876-bd1758c8afd1
# ╠═e5368ab2-eb97-11ea-1755-f5b9bae4e5c3
# ╠═cb4a13ee-eb97-11ea-373e-474e38399b44
# ╟─fe16f312-eb97-11ea-1272-fb0eb517fbe2
# ╠═1059f184-eb98-11ea-1108-a50008b9be0c
# ╠═eb567d8c-122d-11eb-0f13-f7ea2dc1f805
# ╟─37924058-eb98-11ea-3a78-390939122024
# ╠═61681d9c-eb98-11ea-1411-df9bdf7296ae
# ╟─686ef54a-eb98-11ea-1b68-d33aa9873780
# ╠═5b284ec8-eb97-11ea-142f-85ff5053bd2b
# ╟─b0943d46-eb98-11ea-1a87-8d8b461219a8
# ╠═d7308a84-eb98-11ea-2299-95175b5cf27d
# ╟─e17084e2-eb98-11ea-34ef-b1ebfab71041
