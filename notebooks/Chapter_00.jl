### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ e6dd9f40-a007-439f-bafd-0127f2d7b742
using Pkg

# ╔═╡ f5359bf0-28ff-4ff5-8e70-a2bcc687dfc0
Pkg.activate(expanduser("~/.julia/dev/SR2TuringPluto"))

# ╔═╡ bd37f4aa-5e4e-4823-b758-16c7ee2b2fa0
begin
	using CSV, DataFrames
	using GLM
	using StatisticalRethinking
	using StatisticalRethinkingPlots
	using PlutoUI
end

# ╔═╡ 397f31b2-eba9-11ea-3fd5-1704e3a50a96
md"## Chapter_00.jl"

# ╔═╡ eb5b5d50-c274-4d42-a6a3-5c03c7ff1c14
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

# ╔═╡ 2c7baec4-0f38-4050-be86-b54de4112dc7
md" ### Code 0.0"

# ╔═╡ 1b540d64-dc1f-11ea-321a-f1e0d2bab3e8
md"### Code 0.1"

# ╔═╡ 624b8daa-dc1f-11ea-114f-19896cd76aa5
"All models are wrong, but some are useful."

# ╔═╡ 8fa79b3e-dc1f-11ea-1437-d7bfe6cdea51
md"### Code 0.2"

# ╔═╡ c4321a66-de4f-11ea-3449-774edd51eb10
@bind N Slider(1:5, default=3)

# ╔═╡ 35330a6a-dc20-11ea-2929-9b904165e207
md"##### Variable x initially is a StepRange, not a vector. The log.(x) notation `broadcast` the log function to all steprange elements in x and returns a vector."

# ╔═╡ df73cd1c-dc1f-11ea-0e14-c7e381f49c4a
begin
	x = 1:N
	x = x*10
	x = log.(x)
	x = sum(x)
	x = exp(x)
	x = x*10
	x = log(x)
	x = sum(x)
	x = exp(x)
end

# ╔═╡ 975e15f4-122c-11eb-2a4d-ed836c7bade0
md"##### Notice that in Pluto notebooks variables can only be defined once. This is needed in Pluto to make the `reactivity` work (e.g. when updating above slider setting). Bracketing above sequence of assignments with a `begin` and `end` is allowed."

# ╔═╡ b1d6b4ac-dc20-11ea-3237-a515964f9176
md"### Code 0.3"

# ╔═╡ 88e3bb64-dc1e-11ea-1136-bd822afb72e4
[log(0.01^200) 200 * log(0.01)]

# ╔═╡ dc515b15-a11a-4372-9fad-c00d9c508dd5
md"### Code 0.4"

# ╔═╡ 3a3e7ff9-a76c-4766-9dc8-4de5a4a14caa
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame; delim=';')
	df = filter(row -> row[:age] >= 18, df);
end;

# ╔═╡ 5d3acee5-3515-44e7-80a5-a4927f854554
describe(df)

# ╔═╡ b8339211-7cf8-4dd0-8009-05d410373922
md"##### Fit a linear regression of distance on speed."

# ╔═╡ 32f1748c-e8fc-4afc-8d03-d79afa01ad36
m = lm(@formula(height ~ weight), df)

# ╔═╡ c186e7ab-0c18-4a5c-a61f-a061389e21c2
md"##### Estimated coefficients from the model."

# ╔═╡ 4655ab6c-0951-4a3b-bfc8-e991aaf0eff7
coef(m)

# ╔═╡ 88eaa866-b87e-419b-a910-309c07cbde3c
md"##### Plot residuals against height."

# ╔═╡ 9641b455-c5dd-4cd1-9d9c-aebe91139f5f
scatter( df.height, residuals(m), xlab="Height",
  ylab="Model residual values", lab="Model residuals", leg=:bottomright)

# ╔═╡ bc3de68d-e669-4233-96df-9ff82865ccfb
md"### Code 0.5"

# ╔═╡ 14c96b8b-dca7-4d6f-a76a-e9947bf5c3bd
md" ##### See Code 0.0"

# ╔═╡ c4403280-dc20-11ea-0550-5f2f606a9fea
md"## End of Chapter_00.jl"

# ╔═╡ Cell order:
# ╟─397f31b2-eba9-11ea-3fd5-1704e3a50a96
# ╠═eb5b5d50-c274-4d42-a6a3-5c03c7ff1c14
# ╟─2c7baec4-0f38-4050-be86-b54de4112dc7
# ╠═e6dd9f40-a007-439f-bafd-0127f2d7b742
# ╠═f5359bf0-28ff-4ff5-8e70-a2bcc687dfc0
# ╠═bd37f4aa-5e4e-4823-b758-16c7ee2b2fa0
# ╟─1b540d64-dc1f-11ea-321a-f1e0d2bab3e8
# ╠═624b8daa-dc1f-11ea-114f-19896cd76aa5
# ╟─8fa79b3e-dc1f-11ea-1437-d7bfe6cdea51
# ╠═c4321a66-de4f-11ea-3449-774edd51eb10
# ╟─35330a6a-dc20-11ea-2929-9b904165e207
# ╠═df73cd1c-dc1f-11ea-0e14-c7e381f49c4a
# ╟─975e15f4-122c-11eb-2a4d-ed836c7bade0
# ╟─b1d6b4ac-dc20-11ea-3237-a515964f9176
# ╠═88e3bb64-dc1e-11ea-1136-bd822afb72e4
# ╟─dc515b15-a11a-4372-9fad-c00d9c508dd5
# ╠═3a3e7ff9-a76c-4766-9dc8-4de5a4a14caa
# ╠═5d3acee5-3515-44e7-80a5-a4927f854554
# ╟─b8339211-7cf8-4dd0-8009-05d410373922
# ╠═32f1748c-e8fc-4afc-8d03-d79afa01ad36
# ╟─c186e7ab-0c18-4a5c-a61f-a061389e21c2
# ╠═4655ab6c-0951-4a3b-bfc8-e991aaf0eff7
# ╟─88eaa866-b87e-419b-a910-309c07cbde3c
# ╠═9641b455-c5dd-4cd1-9d9c-aebe91139f5f
# ╟─bc3de68d-e669-4233-96df-9ff82865ccfb
# ╟─14c96b8b-dca7-4d6f-a76a-e9947bf5c3bd
# ╟─c4403280-dc20-11ea-0550-5f2f606a9fea
