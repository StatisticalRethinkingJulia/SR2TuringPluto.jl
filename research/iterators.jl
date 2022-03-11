struct Squares
    count::Int
end

Base.iterate(S::Squares, state=1) =
    state > S.count ? nothing : (state*state, state+1)

for item in Squares(9) # or "for item = Squares(9)"
    # body
    display(item)
end

println(25 in Squares(10))

Base.eltype(::Type{Squares}) = Int # Note that this is defined for the type
Base.length(S::Squares) = S.count
Base.sum(S::Squares) = (n = S.count; return n*(n+1)*(2n+1)÷6)

collect(Squares(4)) |> display

