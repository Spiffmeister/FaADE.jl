

struct grid
    domain  :: NamedTuple
    Δ       :: Tuple
    n       :: Tuple

    function grid(𝒟::Array,npts::Int)
        Δ = collect(range(𝒟[1],𝒟[2],length=npts))

        new(𝒟,Δ,npts)
    end
end



