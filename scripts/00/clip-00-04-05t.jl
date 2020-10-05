
using Markdown
using InteractiveUtils

using DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

md"## Clip-00-04-05t.jl"

md"### snippet 0.5"

md"##### Snippet 0.5 is replaced by below lines."

md"### snippet 0.4"

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame; delim=';')
	howell1 = filter(row -> row[:age] >= 18, df);
	first(howell1, 5)
end

md"##### Fit a linear regression of distance on speed."

m = lm(@formula(height ~ weight), howell1)

md"##### Estimated coefficients from the model."

coef(m)

md"##### Plot residuals against height."

scatter( howell1.height, residuals(m), xlab="Height",
  ylab="Model residual values", lab="Model residuals", leg=:bottomright)

md"## End of clip-00-04-05t.jl"

