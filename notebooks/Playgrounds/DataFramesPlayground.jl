### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 58c3936e-c10a-40bd-b8f9-a15804c27e59
using DataFrames

# ╔═╡ 636f1f75-2055-4cd3-95bc-3e4f1885f79c
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


# ╔═╡ 2a70a1c7-054f-43d1-84ac-0c1b737233cb
df = DataFrame(a=1:3, b=4:6)

# ╔═╡ 8993cf19-ae9b-4add-95bc-78207dd6ec9c
md" #### `getindex` (used in df on the right hand side of =)."

# ╔═╡ d2959067-58d9-487f-af6f-2ce2980d63ae
let
	copyOfa = df[:, :a]
	copyOfa .= 0
	copyOfa
end

# ╔═╡ 2434c37c-a7d4-4e1a-985d-7824885fbab4
df

# ╔═╡ f4ef304f-2ca5-4e95-86b0-d83865d74aeb
md" ##### Note that column `:a` has not changed."

# ╔═╡ fb181af4-b916-4aac-b997-1b89dbb99fb5
md" #### `setindex` (used in df on left hand side of =)."

# ╔═╡ 2b49dcd4-4054-4300-b75a-c109b40fc815
begin
	df[:, :a] .= 0
	df
end

# ╔═╡ 957207cf-47ce-4f2d-ba3b-53b609401c8d
md" ##### Now df has changed!"

# ╔═╡ a373977d-c27b-4478-8beb-58b15b7995a9
begin
	df .= 0
	df
end

# ╔═╡ 6111833d-6c40-4de1-9692-93a276a5730a
md" ##### In `getindex` context a user can ask for the following behaviors:

1. get a vector “as is”, without copying; this is achieved with df[!, :a] and df.a;
2. get a vector with copying; that is achieved with df[:, :a];
3. get a view of the vector; this is achieved with view(df, :, :a) or view(df, !, :a)
"

# ╔═╡ 331bde26-9892-4045-8355-85f666514487
md""" ##### In `setindex!` + broadcasting assignment context the user might ask for the following:

1. set an existing vector in-place, this is done by: df[:, :a] = vec and df[:, :a] .= something

2. replace an existing vector without copying: df[!, :a] = vec or df.a = vec

3. replace an existing vector with copying: df[!, :a] .= vec and df[!, :a] .= something

4. add a new vector without copying: df[!, :new] = vec or df.new = vec

5. add a new vector with copying: df[:, :new] = vec or df[:, :new] .= something

or 

6. df[!, :new] .= something

(this last rule breaks things a bit as ! and : behave here the in same way, but it was added for user convenience)
"""

# ╔═╡ 9e510dbd-2e9c-40b3-a5ca-4ab4aa02d019
begin
	a = [7; 8; 9]
	df[:, :a] = a
	df
end

# ╔═╡ 2d15edc9-1f27-4d89-8824-7f6f27385b67
let
	b = [-1; -2; -3]
	df[!, :a] = b
	df
end

# ╔═╡ 5e631910-7a4f-4651-b605-c12966468757
begin
	b = 10:12
	df
end

# ╔═╡ 2b33122b-4e96-42b4-a2f3-83b6f231fd1c
df

# ╔═╡ 9d8a3a01-7c1c-4505-b7d1-42849d44ce3d


# ╔═╡ Cell order:
# ╠═636f1f75-2055-4cd3-95bc-3e4f1885f79c
# ╠═58c3936e-c10a-40bd-b8f9-a15804c27e59
# ╠═2a70a1c7-054f-43d1-84ac-0c1b737233cb
# ╟─8993cf19-ae9b-4add-95bc-78207dd6ec9c
# ╠═d2959067-58d9-487f-af6f-2ce2980d63ae
# ╠═2434c37c-a7d4-4e1a-985d-7824885fbab4
# ╟─f4ef304f-2ca5-4e95-86b0-d83865d74aeb
# ╟─fb181af4-b916-4aac-b997-1b89dbb99fb5
# ╠═2b49dcd4-4054-4300-b75a-c109b40fc815
# ╟─957207cf-47ce-4f2d-ba3b-53b609401c8d
# ╠═a373977d-c27b-4478-8beb-58b15b7995a9
# ╟─6111833d-6c40-4de1-9692-93a276a5730a
# ╟─331bde26-9892-4045-8355-85f666514487
# ╠═9e510dbd-2e9c-40b3-a5ca-4ab4aa02d019
# ╠═2d15edc9-1f27-4d89-8824-7f6f27385b67
# ╠═5e631910-7a4f-4651-b605-c12966468757
# ╠═2b33122b-4e96-42b4-a2f3-83b6f231fd1c
# ╠═9d8a3a01-7c1c-4505-b7d1-42849d44ce3d
