# Clip-03-01t.jl

using DrWatson
@quickactivate "StatisticalRethinkingTuring"

# snippet 3.1

Pr_Positive_Vampire = 0.95
Pr_Positive_Mortal = 0.01
Pr_Vampire = 0.001
Pr_Positive = Pr_Positive_Vampire * Pr_Vampire + Pr_Positive_Mortal * (1 - Pr_Vampire)
Pr_Vampire_Positive = Pr_Positive_Vampire * Pr_Vampire / Pr_Positive
Pr_Vampire_Positive

# End of clip-03-01t.jl
