
"""
    SAT_Robin

Simulatenous approximation term for Robin boundary conditions
`NodeType` is either `Left` or `Right`

TODO: Testing
"""
struct SAT_Robin{T} <: SimultanousApproximationTerm{:Robin}
    type    :: BoundaryConditionType
    side    :: NodeType
    axis    :: Int
    order   :: Int
    ED₁     :: Vector{T}
    RHS     :: Function
    Δx      :: T
    τ       :: T
    a       :: T
    b       :: T

    function SAT_Robin(RHS::Function,a,b,Δx::T,side::NodeType,axis::Int,order::Int) where T
        check_boundary(side)

        τ = SATpenalties(Robin,a,Δx)

        new{T}(Robin,side,axis)
    end
end





function SAT_Robin end
function SAT_Robin(::NodeType{:Left},u::Vector{Float64},Δx::Float64;
        order=2,a=1.0,b=1.0,forcing=false)
    # Get penalties
    τ = SATpenalties(Robin,a,Δx,order)
    SAT = zeros(Float64,order)

    if !forcing
        SAT[1] = b*u[1] + a*BD₁(u,Left,Δx,order)
        return SAT
    else
        SAT[1] = -τ*u[1]
        return SAT
    end
end
function SAT_Robin(::NodeType{:Right},u::Vector{Float64},Δx::Float64;
        order=2,a=1.0,b=1.0,forcing=false)
    # Get penalties
    τ = SATpenalties(Robin,a,Δx,order)
    SAT = zeros(Float64,order)

    if !forcing
        SAT[end] = b*u[end] - a*BD₁(u,Right,Δx,order)
        return SAT
    else
        SAT[end] = τ*u[end]
        return SAT
    end
end




function SAT_Robin(::NodeType{:Left},u::AbstractVector,Δx::Float64,RHS;
        order=2,a=1.0,b=1.0,forcing=false)
    # Get penalties
    τ = SATpenalties(Robin,a,Δx,order)
    SAT = zeros(Float64,order)

    SAT[1] = τ * (b*u[1] + a*BD₁(u,Left,Δx,order) - RHS[1])
    return SAT
end
function SAT_Robin(::NodeType{:Right},u::AbstractVector,Δx::Float64,RHS;
        order=2,a=1.0,b=1.0,forcing=false)
    # Get penalties
    τ = SATpenalties(Robin,a,Δx,order)
    SAT = zeros(Float64,order)

    SAT[end] = τ * (b*u[end] - a*BD₁(u,Right,Δx,order) + RHS[end])
    return SAT
end
