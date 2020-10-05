
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
  	using Turing
	using StatisticalRethinking
	using PlutoUI
end

md"## Fig 2.5t"

md"### This clip is only intended to generate Fig 2.5."

@model globe_toss(n, k) = begin
  theta ~ Beta(1, 1) # prior
  k ~ Binomial(n, theta) # model
  return k, theta
end;

md"### 1. Create a SampleModel object:"

md"##### n will go from 1:9"

@bind n Slider(1:18, default=9)

begin
	k = [1,0,1,1,1,0,1,0,1]       # Sequence observed
	x = range(0, stop=9, length=10)
end

begin
	p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 9)
	dens = Vector{DataFrame}(undef, 10)
	for n in 1:9
		p[n] = plot(xlims=(0.0, 1.0), ylims=(0.0, 3.0), leg=false)
		k = sum([1,0,1,1,1,0,1,0,1,1,0,1,1,1,0,1,0,1][1:n])
		chns = sample(globe_toss(n, k), NUTS(0.65), 1000)
		dfs=DataFrame(chns)
		if n == 1
			hline!([1.0], line=(:dash))
		else
			density!(dens[n][:, :theta], line=(:dash))
		end
		density!(dfs[:, :theta])
		dens[n+1] = dfs

	end
end


plot(p..., layout=(3, 3))

md"## End of Fig2.5t.jl"

