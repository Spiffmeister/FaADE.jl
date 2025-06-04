


"""
    cosinespace(R::Vector{TT},θ::TT,ζ::TT,m::Vector,n::Vector) where TT
Cosine term for cylindrical coordinates
"""
function cosinespace(R::Array{TT},θ::TT,ζ::TT,m::Vector,n::Vector) where TT
    local r = 0.0 :: TT
    # for j in eachindex(m)
        for i in eachindex(n)
            r += R[i]*cos(convert(TT,m[i])*θ - convert(TT,n[i])*ζ)
        end
    # end
    return r
end
function cosinespace(R::Array{TT},θ::AbstractVector{TT},ζ,m,n) where TT
    r = zeros(TT,length(θ))
    for i in eachindex(θ)
        r[i] = cosinespace(R,θ[i],ζ,m,n)
    end
    return r
end

"""
    sinusespace(Z::Vector{TT},θ::TT,ζ::TT,m::Vector,n::Vector) where TT
Sin term for cylindrical coordinates
"""
function sinusespace(Z::Array{TT},θ::TT,ζ::TT,m::Vector{Int},n::Vector{Int}) where TT
    local z = 0.0 :: TT
    # for j in eachindex(m)
        for i in eachindex(n)
            z += Z[i]*sin(m[i]*θ - n[i]*ζ)
        end
    # end
    return z
end
function sinusespace(Z::Array{TT},θ::AbstractVector{TT},ζ,m,n) where TT
    z = zeros(TT,length(θ))
    for i in eachindex(θ)
        z[i] = sinusespace(Z,θ[i],ζ,m,n)
    end
    return z
end




"""
    Torus{TT}
Represents a torus in boundary using Fourier series in `R` and `Z` coordinates as,

``R(\\theta,\\zeta) = \\sum_{i,j} R_{i,j} cos(m_j\\theta - n_i\\zeta), \\qquad Z(\\theta,\\zeta) = \\sum_{i,j} Z_{i,j} sin(m_j\\theta - n_i\\zeta)``

where `Rᵢⱼ` and `Zᵢⱼ` are the Fourier coefficients and `mⱼ` and `nᵢ` are the Fourier modes.

Example:
```julia
Rin = [3e-1]; Zin=[-3e-1]
Rout = [6e-1]; Zout=[-6e-1]

inner = Torus(Rout,Zout,[1],[0])
outer = Torus(Rout,Zout,[1],[0])

meshgrid(inner,outer,0.0,11,21)
```

Coordinates at a given ``\\theta\\in[0,R]`` and ``\\zeta\\in[0,2\\pi)`` can be computed using the call syntax
```julia
inner(θ,ζ)
```
"""
struct Torus{TT}
    R::Array{TT}
    Z::Array{TT}
    m::Vector{Int}
    n::Vector{Int}
end
"""
    (T::Torus)(θ::TT,ζ::TT) where TT
Compute the coordinates of the torus at the given `θ` and `ζ` values.
"""
function (T::Torus)(θ::TT,ζ::TT) where TT
    return [cosinespace(T.R,θ,ζ,T.m,T.n), sinusespace(T.Z,θ,ζ,T.m,T.n)]
end

