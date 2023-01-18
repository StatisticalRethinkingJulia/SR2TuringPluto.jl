### A Pluto.jl notebook ###
# v0.19.3

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

using Markdown
using InteractiveUtils

# ╔═╡ 80652cb5-a82f-4a28-8212-4bb42edc2cfe
using Pkg, DrWatson

# ╔═╡ a76198d0-28ba-4aaa-83d3-e5417f7bcd66
begin
	using DataFrames
	using StatsBase
	using StatsPlots
	using Random
	using Distributions
end

# ╔═╡ a532c481-bcca-4c6f-b297-a2e18e1cebe1
default(label=false);

# ╔═╡ f61eb194-60e0-423d-9cf9-7200b3985ce3
md" ## 10.1 Maximum entropy."

# ╔═╡ 8cf0ad12-c812-4b9f-8091-c24ea45944bf
md" #### Code 10.1"

# ╔═╡ 2322642e-c12c-49ee-bcff-abc0b7b06e6f
p = DataFrame(
    :A => [0, 0, 10, 0, 0],
    :B => [0, 1, 8, 1, 0],
    :C => [0, 2, 6, 2, 0],
    :D => [1, 2, 4, 2, 1],
    :E => [2, 2, 2, 2, 2],
)

# ╔═╡ 3fa5f5b5-8942-469d-b1f1-b8fa128e7ccc
md" #### Code 10.2"

# ╔═╡ a7a59f9a-8f06-46fd-9142-26f8737331b9
p_norm = mapcols(c -> c ./ sum(c), p)

# ╔═╡ 02ef21ff-bd5e-46c2-8d8e-a14fd7e3fcd4
md" #### Code 10.3"

# ╔═╡ 02e6f217-bb5e-4fac-b870-a46741610fce
ent_vals = mapcols(entropy, p_norm)

# ╔═╡ b26bc5ae-81e9-4cb7-a051-402d5371cc52
md" #### Code 10.4"

# ╔═╡ 5a81f408-ef16-420a-aff2-29c411c970bf
begin
	ways = [1, 90, 1260, 37800, 113400]
	logwayspp = log.(ways)/10
	
	txt = text.(names(ent_vals), :bottom, :right, 11)
	scatter(logwayspp, collect(ent_vals[1,:]), txt=txt, xlab="log(ways) per pebble", ylab="entropy")
end

# ╔═╡ b36f5641-ac74-4a85-b7a9-c2230c2cb1d2
md" #### Code 10.5"

# ╔═╡ 215fe4f3-2163-41f1-b716-470e0934bf6b
p2 = [
    [1/4, 1/4, 1/4, 1/4],
    [2/6, 1/6, 1/6, 2/6],
    [1/6, 2/6, 2/6, 1/6],
    [1/8, 4/8, 2/8, 1/8],
]

# ╔═╡ 4cbca251-6a44-4945-b6de-ae2f6e9fe619
map(x -> sum(x .* [0, 1, 1, 2]), p2)

# ╔═╡ fe87a6c6-a5a5-4e77-bc3c-72ffe44e3d32
md" #### Code 10.6"

# ╔═╡ 023164a1-d894-4e9c-8ba5-34a73133e12b
# Could be simplified with just `map(entropy, p)`
# compute entropy of each distribution
map(x -> -sum(x .* log.(x)), p2)

# ╔═╡ d743025d-ad53-49a1-9998-4518a9e43cd7
md" #### Code 10.7"

# ╔═╡ f0e0b700-4b5c-444d-9cc3-6bbb5dd0c84b
begin
	p3 = 0.7
	A = [
	    (1-p3)^2, p3*(1-p3), (1-p3)*p3, p3^2,
	]
end

# ╔═╡ a29975c2-499b-44f6-a849-e87416cca9c2
md" #### Code 10.8"

# ╔═╡ a28c50ca-2c98-4db2-aa98-a0f51b0e57ad
-sum(A.*log.(A))

# ╔═╡ 7ad2bc63-9933-4e86-b334-96307142ac42
md" #### Code 10.9"

# ╔═╡ a87675f7-533e-4664-84a9-323e14309296
function sim_p(G::Float64 = 1.4)
    p = rand(Uniform(), 3)
    x4 = (G * sum(p) - p[2] - p[3])/(2-G)
    push!(p, x4)
    p ./= sum(p)
    (entropy(p), p)
end

# ╔═╡ 96db954a-c345-4451-bd3e-40ad6fdd2d7d
md" #### Code 10.10"

# ╔═╡ 198b1ae7-2b03-4f8c-900e-855525762a34
begin
	Random.seed!(1)
	cnt = 10^5
	H = [sim_p(1.4) for _ ∈ 1:cnt];
	density(first.(H))
end

# ╔═╡ 68f78d0f-2471-4fa7-96d8-f02679adda0f
md" #### Code 10.11"

# ╔═╡ fefb3fc2-3b0f-4344-b590-53cba794faa0
begin
	entropies = first.(H)
	distributions = last.(H);
end;

# ╔═╡ d7d1cf81-4903-4924-86aa-1b10e87d29b4
md" #### Code 10.12"

# ╔═╡ 5773f70f-0d01-46c6-934d-90d39e0a792e
maximum(entropies)

# ╔═╡ fb55b28a-8930-4b00-b1ca-e8007662ca7c
md" #### Code 10.13"

# ╔═╡ 342961de-a5e0-42d1-b127-f9d825c2b265
distributions[findmax(entropies)[2]]

# ╔═╡ b29892a0-9f11-486a-a708-4a565c318a0d
md" ## 10.2 Generalized linear models."

# ╔═╡ 5b7f5df0-05aa-4bcc-9a5a-911694ce4f85
# No code here, wheee!

# ╔═╡ Cell order:
# ╠═80652cb5-a82f-4a28-8212-4bb42edc2cfe
# ╠═a76198d0-28ba-4aaa-83d3-e5417f7bcd66
# ╠═a532c481-bcca-4c6f-b297-a2e18e1cebe1
# ╟─f61eb194-60e0-423d-9cf9-7200b3985ce3
# ╟─8cf0ad12-c812-4b9f-8091-c24ea45944bf
# ╠═2322642e-c12c-49ee-bcff-abc0b7b06e6f
# ╟─3fa5f5b5-8942-469d-b1f1-b8fa128e7ccc
# ╠═a7a59f9a-8f06-46fd-9142-26f8737331b9
# ╠═02ef21ff-bd5e-46c2-8d8e-a14fd7e3fcd4
# ╠═02e6f217-bb5e-4fac-b870-a46741610fce
# ╠═b26bc5ae-81e9-4cb7-a051-402d5371cc52
# ╠═5a81f408-ef16-420a-aff2-29c411c970bf
# ╠═b36f5641-ac74-4a85-b7a9-c2230c2cb1d2
# ╠═215fe4f3-2163-41f1-b716-470e0934bf6b
# ╠═4cbca251-6a44-4945-b6de-ae2f6e9fe619
# ╟─fe87a6c6-a5a5-4e77-bc3c-72ffe44e3d32
# ╠═023164a1-d894-4e9c-8ba5-34a73133e12b
# ╟─d743025d-ad53-49a1-9998-4518a9e43cd7
# ╠═f0e0b700-4b5c-444d-9cc3-6bbb5dd0c84b
# ╟─a29975c2-499b-44f6-a849-e87416cca9c2
# ╠═a28c50ca-2c98-4db2-aa98-a0f51b0e57ad
# ╟─7ad2bc63-9933-4e86-b334-96307142ac42
# ╠═a87675f7-533e-4664-84a9-323e14309296
# ╟─96db954a-c345-4451-bd3e-40ad6fdd2d7d
# ╠═198b1ae7-2b03-4f8c-900e-855525762a34
# ╟─68f78d0f-2471-4fa7-96d8-f02679adda0f
# ╠═fefb3fc2-3b0f-4344-b590-53cba794faa0
# ╟─d7d1cf81-4903-4924-86aa-1b10e87d29b4
# ╠═5773f70f-0d01-46c6-934d-90d39e0a792e
# ╟─fb55b28a-8930-4b00-b1ca-e8007662ca7c
# ╠═342961de-a5e0-42d1-b127-f9d825c2b265
# ╟─b29892a0-9f11-486a-a708-4a565c318a0d
# ╠═5b7f5df0-05aa-4bcc-9a5a-911694ce4f85
