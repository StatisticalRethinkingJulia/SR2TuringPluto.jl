### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ d388f90e-7530-415a-99f1-32472a77db9c
using Pkg

# ╔═╡ b7c94aeb-b8ca-49bd-8ece-dbbf6611e086
begin
	using Distributions
	using StatsPlots
	using StatsBase
	using KernelDensity
	using StatisticalRethinking
end

# ╔═╡ 81f65c1f-c44b-4e71-97e5-ae62a38c9b55
default(labels=false)

# ╔═╡ 709a8185-d0d7-4367-866f-acf0aa00abe4
md"### Code 3.1"

# ╔═╡ d6aed17b-5b4c-4ed3-86c9-d72e3d3f2d57
let
	Pr_Positive_Vampire = 0.95
	Pr_Positive_Mortal = 0.01
	Pr_Vampire = 0.001
	tmp = Pr_Positive_Vampire * Pr_Vampire
	Pr_Positive = tmp + Pr_Positive_Mortal * (1 - Pr_Vampire)
	Pr_Vampire_Positive = tmp / Pr_Positive
	Pr_Vampire_Positive
end

# ╔═╡ 10429558-28dc-477a-80cc-b0cf31976fac
md"### Code 3.2 & 3.3"

# ╔═╡ 4da7c481-ce34-475e-af71-007d08cfcb6f
begin
	size = 1000
	p_grid = range(0, 1; length=size)
	prob_p = repeat([1.0], size);
	prob_data = [pdf(Binomial(9, p), 6) for p in p_grid];
	posterior = prob_data .* prob_p
	posterior /= sum(posterior);

	samples_count = 10_000
	cat = Categorical(posterior);
	indices = rand(cat, samples_count)
	samples = p_grid[indices];
end

# ╔═╡ 83cc1cc6-5989-41ba-8680-12b03c540288
md"### Code 3.4"

# ╔═╡ 79859429-7826-480d-8ef0-47e77257c644
scatter(samples; alpha=0.2)

# ╔═╡ 6ae0be1b-c758-4fef-abc8-f23026eb3a87
md"### Code 3.5"

# ╔═╡ d593734f-6a5d-4f20-bf52-c0f46f83253b
density(samples)

# ╔═╡ 55066707-fa50-471c-a624-1fd40723687d
md"### Code 3.6"

# ╔═╡ b013128d-750c-48e5-b11c-87e7faa28022
sum(posterior[p_grid .< 0.5])

# ╔═╡ 984f45be-5c69-425f-be06-b34dc2955b68
md"### Code 3.7"

# ╔═╡ 1fe37c5c-c2fc-4198-8f9c-0e51d51006bc
sum(samples .< 0.5) / samples_count

# ╔═╡ 231460fd-f576-45a8-a587-7d285439ff02
md"### Code 3.8"

# ╔═╡ 3e528812-93ab-4ffc-bb56-fa40f9cefb8c
sum(@. (samples > 0.5) & (samples < 0.75)) / samples_count

# ╔═╡ 360b2cb7-ab6a-4f17-851f-1804ebd6af4f
md"### Code 3.9"

# ╔═╡ 53da878f-3bde-45ed-a2da-ee0dc03912a0
quantile(samples, 0.8)

# ╔═╡ 322d5b46-b204-42e6-bf4d-bc8a8aa9e091
md"### Code 3.10"

# ╔═╡ 23fe9c96-fa08-47b5-9ce6-97fdcd5915ed
quantile(samples, [0.1, 0.9])

# ╔═╡ ccea74d5-f3b5-4b5d-98f7-a42623ed789d
md"### Code 3.11"

# ╔═╡ 58637734-3529-4544-85ef-f982a17b5894
begin
	prob_data2 = [pdf(Binomial(3, p), 3) for p in p_grid];
	posterior2 = prob_data2 .* prob_p
	posterior2 /= sum(posterior2)

	cat2 = Categorical(posterior2);
	samples2 = p_grid[rand(cat2, samples_count)];
end

# ╔═╡ e0d0764e-4d9d-4009-be03-e882c8d782d3
md"### Code 3.12"

# ╔═╡ 90986735-3d78-40d6-8ec9-217fea7bd8df
percentile(samples2, [25, 75])

# ╔═╡ 1037f18a-1aae-4532-b9b1-b2c3dddd8d89
md"### Code 3.13"

# ╔═╡ 9540fe23-6e1d-4b88-ac41-10f523286e0a
hpdi(samples2, alpha=0.5)

# ╔═╡ 84b0ee54-4b24-4283-9798-7ec63552f626
md"### Code 3.14"

# ╔═╡ c47a769d-fb74-4397-bd49-594a6281185f
p_grid[argmax(posterior2)]

# ╔═╡ 0c379759-9e6d-4e78-a189-4b0907906e12
md"### Code 3.15"

# ╔═╡ 969ad003-e132-45c9-a7a8-3bbedbc50f43
k = kde(samples2, bandwidth=0.01)

# ╔═╡ ae92a07a-9473-4ac7-8e5b-29023e9054a1
k.x[argmax(k.density)]

# ╔═╡ ea838454-aa3f-4d26-9c59-871c3779580f
md"### Code 3.16"

# ╔═╡ 2643b2a8-41fa-4a8c-b154-70c253b5c414
mean(samples2), median(samples2)

# ╔═╡ 9a1863c0-8c7e-4962-882c-6ddb4968a4d6
md"### Code 3.17"

# ╔═╡ e938bb1c-d47a-4627-a850-0ac255892482
sum(@. posterior2 * abs(0.5 - p_grid))

# ╔═╡ 37a8aeac-1942-4b95-8ead-98d86c03f2bc
md"### Code 3.18"

# ╔═╡ 678b1bf3-f137-46a7-a558-1a8c22964f89
loss = map(d -> sum(@. posterior2 * abs(d - p_grid)), p_grid);

# ╔═╡ 135eaf3d-e370-4f0b-99ba-55c30884d23a
md"### Code 3.19"

# ╔═╡ 2f3c79b4-0b8e-414a-9287-0c74b95e176c
p_grid[argmin(loss)]

# ╔═╡ 72a788c7-98bb-4b76-becd-f230c029266f
md"### Code 3.20"

# ╔═╡ 43d5afea-726c-4948-a66b-6c7660e3d877
[pdf(Binomial(2, 0.7), n) for n ∈ 0:2]

# ╔═╡ ce3b7ea6-1e70-4fe4-833f-8d730c3d10f5
md"### Code 3.21"

# ╔═╡ e7092c37-92e8-403a-ba4e-f109ac90cfbe
rand(Binomial(2, 0.7))

# ╔═╡ 1be7df70-fd63-49cf-90dd-3ecb14a36f79
md"### Code 3.22"

# ╔═╡ 599e3517-7b92-4fb3-ad8d-85105e1ccdf5
s = rand(Binomial(2, 0.7), 10)

# ╔═╡ 2d424f57-18e6-4e72-875c-8e63174bc210
md"### Code 3.23"

# ╔═╡ b7082a08-7300-45d7-9f55-0b151b779be3
let
	dummy_w = rand(Binomial(2, 0.7), 100_000)
	proportions(dummy_w)  # or counts(dummy_w)/100000
end

# ╔═╡ 7ed24a4c-5323-4916-a36b-34b9563c37a4
md"### Code 3.24"

# ╔═╡ 74bf9666-b9ce-4413-bdeb-6427a2901c2c
let
	dummy_w = rand(Binomial(9, 0.7), 100_000)
	histogram(dummy_w; xlabel="dummy water count", ylabel="Frequency")
end

# ╔═╡ dd81f64a-8925-4bdf-9ee1-c4b435db2b0a
md"### Code 3.25"

# ╔═╡ c153f25c-2a07-43bf-8f15-ab02c4793ccc
let
	w = rand(Binomial(9, 0.6), 10_000)
end

# ╔═╡ 61583640-1f75-4a64-b6ad-c220d2994827
md"### Code 3.26"

# ╔═╡ 448335c1-fb2f-4662-a4eb-c4276e0f2156
let
	w = [rand(Binomial(9, p)) for p in samples]
end

# ╔═╡ Cell order:
# ╠═d388f90e-7530-415a-99f1-32472a77db9c
# ╠═b7c94aeb-b8ca-49bd-8ece-dbbf6611e086
# ╠═81f65c1f-c44b-4e71-97e5-ae62a38c9b55
# ╟─709a8185-d0d7-4367-866f-acf0aa00abe4
# ╠═d6aed17b-5b4c-4ed3-86c9-d72e3d3f2d57
# ╟─10429558-28dc-477a-80cc-b0cf31976fac
# ╠═4da7c481-ce34-475e-af71-007d08cfcb6f
# ╟─83cc1cc6-5989-41ba-8680-12b03c540288
# ╠═79859429-7826-480d-8ef0-47e77257c644
# ╟─6ae0be1b-c758-4fef-abc8-f23026eb3a87
# ╠═d593734f-6a5d-4f20-bf52-c0f46f83253b
# ╟─55066707-fa50-471c-a624-1fd40723687d
# ╠═b013128d-750c-48e5-b11c-87e7faa28022
# ╟─984f45be-5c69-425f-be06-b34dc2955b68
# ╠═1fe37c5c-c2fc-4198-8f9c-0e51d51006bc
# ╟─231460fd-f576-45a8-a587-7d285439ff02
# ╠═3e528812-93ab-4ffc-bb56-fa40f9cefb8c
# ╟─360b2cb7-ab6a-4f17-851f-1804ebd6af4f
# ╠═53da878f-3bde-45ed-a2da-ee0dc03912a0
# ╟─322d5b46-b204-42e6-bf4d-bc8a8aa9e091
# ╠═23fe9c96-fa08-47b5-9ce6-97fdcd5915ed
# ╟─ccea74d5-f3b5-4b5d-98f7-a42623ed789d
# ╠═58637734-3529-4544-85ef-f982a17b5894
# ╟─e0d0764e-4d9d-4009-be03-e882c8d782d3
# ╠═90986735-3d78-40d6-8ec9-217fea7bd8df
# ╟─1037f18a-1aae-4532-b9b1-b2c3dddd8d89
# ╠═9540fe23-6e1d-4b88-ac41-10f523286e0a
# ╟─84b0ee54-4b24-4283-9798-7ec63552f626
# ╠═c47a769d-fb74-4397-bd49-594a6281185f
# ╟─0c379759-9e6d-4e78-a189-4b0907906e12
# ╠═969ad003-e132-45c9-a7a8-3bbedbc50f43
# ╠═ae92a07a-9473-4ac7-8e5b-29023e9054a1
# ╟─ea838454-aa3f-4d26-9c59-871c3779580f
# ╠═2643b2a8-41fa-4a8c-b154-70c253b5c414
# ╟─9a1863c0-8c7e-4962-882c-6ddb4968a4d6
# ╠═e938bb1c-d47a-4627-a850-0ac255892482
# ╟─37a8aeac-1942-4b95-8ead-98d86c03f2bc
# ╠═678b1bf3-f137-46a7-a558-1a8c22964f89
# ╟─135eaf3d-e370-4f0b-99ba-55c30884d23a
# ╠═2f3c79b4-0b8e-414a-9287-0c74b95e176c
# ╟─72a788c7-98bb-4b76-becd-f230c029266f
# ╠═43d5afea-726c-4948-a66b-6c7660e3d877
# ╟─ce3b7ea6-1e70-4fe4-833f-8d730c3d10f5
# ╠═e7092c37-92e8-403a-ba4e-f109ac90cfbe
# ╟─1be7df70-fd63-49cf-90dd-3ecb14a36f79
# ╠═599e3517-7b92-4fb3-ad8d-85105e1ccdf5
# ╟─2d424f57-18e6-4e72-875c-8e63174bc210
# ╠═b7082a08-7300-45d7-9f55-0b151b779be3
# ╟─7ed24a4c-5323-4916-a36b-34b9563c37a4
# ╠═74bf9666-b9ce-4413-bdeb-6427a2901c2c
# ╟─dd81f64a-8925-4bdf-9ee1-c4b435db2b0a
# ╠═c153f25c-2a07-43bf-8f15-ab02c4793ccc
# ╟─61583640-1f75-4a64-b6ad-c220d2994827
# ╠═448335c1-fb2f-4662-a4eb-c4276e0f2156
