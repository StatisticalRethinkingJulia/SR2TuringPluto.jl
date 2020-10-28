### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ ddad81d2-184a-11eb-0a24-93759b4a0d17
using Pkg, DrWatson

# ╔═╡ ddd64660-184a-11eb-341e-a1896dcf506b
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing 
	using StatisticalRethinking
	Turing.turnprogress(false);
end

# ╔═╡ f771d004-1849-11eb-308a-fdd834f7e6b1
md"## Predict.jl"

# ╔═╡ ddd6c838-184a-11eb-3285-b704caa74767
begin
	N = 200
	df = DataFrame()
	df.weight = rand(Uniform(20, 90), N)
	f(x) = mean(df.weight) + 1.6x + rand(Normal(0, 10))
	df.height = f.(df.weight)
	Text(precis(df; io=String))
end

# ╔═╡ dddd83d0-184a-11eb-3b52-13e587ea984f
@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* weights
    for i in eachindex(heights)
        heights[i] ~ Normal(μ[i], σ)
    end
end

# ╔═╡ ddeb0d66-184a-11eb-0107-e5532d248439
begin
	m4_3t = m4_3(df.weight, df.height)
	chns4_3t = sample(m4_3t, NUTS(100, 0.65), N);
end

# ╔═╡ ddecc1a6-184a-11eb-2f79-bfe4c56ee921
md"##### Use DataFrame for easy handling."

# ╔═╡ b6ba87d8-184c-11eb-2a47-5f926e636eef
begin
	chns4_3t_df = DataFrame(Array(chns4_3t, :parameters), names(chns4_3t, :parameters))
	Text(precis(chns4_3t_df; io=String))
end

# ╔═╡ ddf925a4-184a-11eb-1556-f3b5938ba663
Particles(chns4_3t[[:a, :b, :σ]])

# ╔═╡ f8b39dc2-1875-11eb-1280-bf16e96a59c5
Particles(chns4_3t_df)

# ╔═╡ de03d602-184a-11eb-093e-a1747b291cc4
md"##### Prediction part."

# ╔═╡ de131f84-184a-11eb-05ec-152d51bde0e2
begin
	x_test = [20, 40, 50, 60, 80]
	m_test = m4_3(x_test, fill(missing, length(x_test)))
	predictions = Turing.Inference.predict(m_test, chns4_3t)
end

# ╔═╡ de14ddf8-184a-11eb-0c1d-b38c9a3d0f49
md"##### Convert to a DataFrame for easy handling."

# ╔═╡ de220e38-184a-11eb-2904-2fd0143ae1fe
begin
	pred_df = DataFrame(predictions)[:, 3:end]
	Text(precis(pred_df; io=String))
end

# ╔═╡ de2b0542-184a-11eb-3bdf-59c61d532a2e
md"#### Get the mean predicted values."

# ╔═╡ de32cf16-184a-11eb-038a-1d7b097ff64f
begin
	pred_summary_df = DataFrame()
	pred_summary_df.parameters = names(predictions)
	pred_summary_df.mean = mean(predictions)[:, 2]
	pred_summary_df.std = vec(std(Array(group(predictions, :heights)); dims = 1))
	pred_summary_df.lower = hpd(predictions; alpha=0.11)[:, 2]
	pred_summary_df.upper = hpd(predictions; alpha=0.11)[:, 3]
	pred_summary_df
end

# ╔═╡ de3a80c4-184a-11eb-221c-5b6199a3f0e1
begin
	scatter(df.weight, df.height, lab="Observations",leg=:bottomright)
	scatter!(x_test, pred_summary_df.mean,
    markersize=5, color=:red, lab="Predictions at x_test")
	x = range(minimum(df.weight), stop=maximum(df.weight), length=N)
	g(x) = mean(chns4_3t[:a]) + mean(chns4_3t[:b]) * x
	plot!(x, g.(x), lab="Regression line", leg=:bottomright)
end

# ╔═╡ de42317c-184a-11eb-0f91-87634b6f5007
md"## End of predict.jl"

# ╔═╡ Cell order:
# ╟─f771d004-1849-11eb-308a-fdd834f7e6b1
# ╠═ddad81d2-184a-11eb-0a24-93759b4a0d17
# ╠═ddd64660-184a-11eb-341e-a1896dcf506b
# ╠═ddd6c838-184a-11eb-3285-b704caa74767
# ╠═dddd83d0-184a-11eb-3b52-13e587ea984f
# ╠═ddeb0d66-184a-11eb-0107-e5532d248439
# ╟─ddecc1a6-184a-11eb-2f79-bfe4c56ee921
# ╠═b6ba87d8-184c-11eb-2a47-5f926e636eef
# ╠═ddf925a4-184a-11eb-1556-f3b5938ba663
# ╠═f8b39dc2-1875-11eb-1280-bf16e96a59c5
# ╟─de03d602-184a-11eb-093e-a1747b291cc4
# ╠═de131f84-184a-11eb-05ec-152d51bde0e2
# ╟─de14ddf8-184a-11eb-0c1d-b38c9a3d0f49
# ╠═de220e38-184a-11eb-2904-2fd0143ae1fe
# ╟─de2b0542-184a-11eb-3bdf-59c61d532a2e
# ╠═de32cf16-184a-11eb-038a-1d7b097ff64f
# ╠═de3a80c4-184a-11eb-221c-5b6199a3f0e1
# ╟─de42317c-184a-11eb-0f91-87634b6f5007
