### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 28f7bd2f-3208-4c61-ad19-63b11dd56d30
using Pkg

# ╔═╡ 2846bc48-7972-49bc-8233-80c7ea3326e6
begin
	using DataFrames
    using RegressionAndOtherStories: reset_selected_notebooks_in_notebooks_df!
end

# ╔═╡ 970efecf-9ae7-4771-bff0-089202b1ff1e
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 0%);
    	padding-right: max(160px, 30%);
	}
</style>
"""

# ╔═╡ d98a3a0a-947e-11ed-13a2-61b5b69b4df5
notebook_files = [
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_00.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_02.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_03.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_04_part_1.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_04_part_2.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_04_part_3.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_05_1_Spurious_Associations.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_05_2_Masked_Relationships.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_05_3_Categorical_Variables.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_06.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_07.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_07_issues.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_08.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_09.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_10.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_11.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_12.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_13.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Chapter_14.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Playgrounds/DataFramesMiniLanguage.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Playgrounds/DataFramesPlayground.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Playgrounds/IteratorPlayground.jl",
    "~/.julia/dev/SR2TuringPluto/notebooks/Playgrounds/TuringGuide.jl",
	"~/.julia/dev/SR2TuringPluto/notebooks/Maintenance/Notebook-to-reset-SR2TuringPluto-jl-notebooks.jl"
];

# ╔═╡ 0f10a758-e442-4cd8-88bc-d82d8de97ede
notebooks_df = DataFrame(
    file = notebook_files,
    reset = repeat([false], length(notebook_files)),
	done = repeat([false], length(notebook_files))
)

# ╔═╡ a4207232-61eb-4da7-8629-1bcc670ab524
notebooks_df.reset .= true;

# ╔═╡ 722d4847-2458-4b23-b6a0-d1c321710a2a
notebooks_df

# ╔═╡ 9d94bebb-fc41-482f-8759-cdf224ec71fb
reset_selected_notebooks_in_notebooks_df!(notebooks_df)

# ╔═╡ 88720478-7f64-4852-8683-6be50793666a
notebooks_df

# ╔═╡ Cell order:
# ╠═28f7bd2f-3208-4c61-ad19-63b11dd56d30
# ╠═2846bc48-7972-49bc-8233-80c7ea3326e6
# ╠═970efecf-9ae7-4771-bff0-089202b1ff1e
# ╠═d98a3a0a-947e-11ed-13a2-61b5b69b4df5
# ╠═0f10a758-e442-4cd8-88bc-d82d8de97ede
# ╠═a4207232-61eb-4da7-8629-1bcc670ab524
# ╠═722d4847-2458-4b23-b6a0-d1c321710a2a
# ╠═9d94bebb-fc41-482f-8759-cdf224ec71fb
# ╠═88720478-7f64-4852-8683-6be50793666a
