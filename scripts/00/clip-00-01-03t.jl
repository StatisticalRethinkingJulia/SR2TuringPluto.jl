
using Markdown
using InteractiveUtils

macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using PlutoUI
end

md"## Clip-00-01-03t.jl"

md"### snippet 0.1"

"All models are wrong, but some are useful."

md"### snippet 0.2"

@bind N Slider(1:5, default=3)

md"##### Variable x initially is a StepRange, not a vector. The log.(x) notation `broadcast` the log function to all steprange elements in x and returns a vector."

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

md"##### Notice that in Pluto notebooks varables can only be defined once. This is needed in Pluto to make the `reactivity` work (e.g. when updating above slider setting). Bracketing above sequence of assignments with a `begin` and `end` is allowed."

md"### snippet 0.3"

[log(0.01^200) 200 * log(0.01)]

md"## End of clip-00-01-03t.jl"

