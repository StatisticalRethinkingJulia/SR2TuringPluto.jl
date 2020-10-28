
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

md"## Broadcasting.jl"

p_grid = range(0, 1, length = 5)

prior = ones(5)

likelihood = pdf.(Binomial.(9, p_grid), 6)

posterior = likelihood .* prior

posterior ./= sum(posterior)

md"##### Broadcasting is another important difference between Julia and R.

You don't have automatic
broadcasting (making functions work the same whether you run it on a single
number or on an array of numbers). However, fixing this is very easy,
you only have to annotate your function call with a dot:
  `f(single_number)` -> `f.(array_of_numbers)`.

That means the translation of the line `likelihood = pdf.(Binomial.(9, p_grid), 6)` goes something like:
1. First or every value in `p_grid` make a binomial distribution with that p value.
2. For every distribution then take the pdf of 6.

The same in the next line, this is elementwise multiplication.

For more infos read [this](https://julialang.org/blog/2017/01/moredots)."

md"## End of broadcasting.jl"

