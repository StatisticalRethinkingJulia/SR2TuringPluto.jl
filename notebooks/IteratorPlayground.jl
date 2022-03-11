### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 4fc964c9-4690-412f-86ca-0ad975e5dc44
struct Squares
    count::Int
end

# ╔═╡ 2a1988c7-e51f-4a8f-b8c0-f71a744b3c68
begin
	Base.iterate(S::Squares, state=1) =
    	state > S.count ? nothing : (state*state, state+1)
			Base.eltype(::Type{Squares}) = Int # Note that this is defined for the type
	Base.length(S::Squares) = S.count
	Base.sum(S::Squares) = (n = S.count; return n*(n+1)*(2n+1)÷6)
	Base.iterate(rS::Iterators.Reverse{Squares}, state=rS.itr.count) = state < 1 ? 		nothing : (state*state, state-1)
	
	function Base.getindex(S::Squares, i::Int)
		1 <= i <= S.count || throw(BoundsError(S, i))
		return i*i
	end
	
	Base.firstindex(S::Squares) = 1
	Base.lastindex(S::Squares) = length(S)
	
	Base.getindex(S::Squares, i::Number) = S[convert(Int, i)]
	Base.getindex(S::Squares, I) = [S[i] for i in I]
end

# ╔═╡ 99e4bdcd-696c-4a87-bd2b-9445b3f32bff
begin
	items = Int[]
	for item in Squares(9) # or "for item = Squares(9)"
	    # body
	    append!(items, item)
	end
	items
end

# ╔═╡ ba0a922b-37fd-4608-a855-0b2fb52e1e3b
25 in Squares(10)

# ╔═╡ fe592470-cca8-4f18-a7c1-ae72ea478067
collect(Squares(8))

# ╔═╡ 13631280-0ee6-466d-b46e-af49d8478702
sum(Squares(3))

# ╔═╡ 542e5321-e4f4-4991-b888-c2b373e07be1
collect(Iterators.reverse(Squares(4)))

# ╔═╡ e0d63466-3c0d-43f1-bc72-09c15e361edd
Squares(100)[23]

# ╔═╡ 41f90116-82f9-4a64-a300-8e9f119cb6f8
Squares(23)[end]

# ╔═╡ 0a09c030-25f0-47ca-ad3b-ad745c1e3052
Squares(10)[[3, 4., 5]]

# ╔═╡ ba97acaf-8b4c-4743-9277-74de9db40ac4
begin
	struct SquaresVector <: AbstractArray{Int, 1}
		count::Int
	end
	Base.size(S::SquaresVector) = (S.count,)
	Base.IndexStyle(::Type{<:SquaresVector}) = IndexLinear()
	Base.getindex(S::SquaresVector, i::Int) = i*i
end

# ╔═╡ fd5337b3-9ef4-4ace-b69f-cfbf036c58d4
s = SquaresVector(4)

# ╔═╡ 320f9fb8-be12-4287-aed7-1847329cc679
s[s .> 8]

# ╔═╡ eb629fdf-7810-4ce4-9133-6825f669ac2b
s + s

# ╔═╡ e50a56db-3eaa-4323-8ae6-b36b716c96e6
sin.(s)

# ╔═╡ a702d978-462c-4879-859d-5a92ba4b6c83
IndexStyle(SquaresVector)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0-DEV"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╠═4fc964c9-4690-412f-86ca-0ad975e5dc44
# ╠═2a1988c7-e51f-4a8f-b8c0-f71a744b3c68
# ╠═99e4bdcd-696c-4a87-bd2b-9445b3f32bff
# ╠═ba0a922b-37fd-4608-a855-0b2fb52e1e3b
# ╠═fe592470-cca8-4f18-a7c1-ae72ea478067
# ╠═13631280-0ee6-466d-b46e-af49d8478702
# ╠═542e5321-e4f4-4991-b888-c2b373e07be1
# ╠═e0d63466-3c0d-43f1-bc72-09c15e361edd
# ╠═41f90116-82f9-4a64-a300-8e9f119cb6f8
# ╠═0a09c030-25f0-47ca-ad3b-ad745c1e3052
# ╠═ba97acaf-8b4c-4743-9277-74de9db40ac4
# ╠═fd5337b3-9ef4-4ace-b69f-cfbf036c58d4
# ╠═320f9fb8-be12-4287-aed7-1847329cc679
# ╠═eb629fdf-7810-4ce4-9133-6825f669ac2b
# ╠═e50a56db-3eaa-4323-8ae6-b36b716c96e6
# ╠═a702d978-462c-4879-859d-5a92ba4b6c83
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
