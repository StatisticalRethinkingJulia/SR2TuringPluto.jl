
using Markdown
using InteractiveUtils

using Pkg, DrWatson

md"## Clip-03-01t.jl"

@quickactivate "StatisticalRethinkingTuring"

md"### snippet 3.1"

begin
	Pr_Positive_Vampire = 0.95
	Pr_Positive_Mortal = 0.01
	Pr_Vampire = 0.001
	Pr_Positive = Pr_Positive_Vampire * Pr_Vampire + Pr_Positive_Mortal * (1 - Pr_Vampire)
	Pr_Vampire_Positive = Pr_Positive_Vampire * Pr_Vampire / Pr_Positive
	Pr_Vampire_Positive
end

md"## End of clip-03-01t.jl"

