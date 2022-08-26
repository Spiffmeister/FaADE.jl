
"""
    SAT_Robin

1. SAT_Robin(::NodeType{:Left},u::Vector{Float64},Δx::Float64;
    order=2,a=1.0,b=1.0,forcing=false)
2. SAT_Robin(::NodeType{:Left},u::AbstractVector,Δx::Float64,RHS;
    order=2,a=1.0,b=1.0,forcing=false)


Simulatenous approximation term for Robin boundary conditions
`NodeType` is either `Left` or `Right`
"""
function SAT_Robin end
function SAT_Robin(::NodeType{:Left},u::Vector{Float64},Δx::Float64;
        order=2,a=1.0,b=1.0,forcing=false)
    # Get penalties
    τ = SATpenalties(Robin,a,Δx,order)
    SAT = zeros(Float64,order)

    if !forcing
        SAT[1] = b*u[1] + a*BDₓ(u,Left,Δx,order)
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
        SAT[end] = b*u[end] - a*BDₓ(u,Right,Δx,order)
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

    SAT[1] = τ * (b*u[1] + a*BDₓ(u,Left,Δx,order) - RHS[1])
    return SAT
end
function SAT_Robin(::NodeType{:Right},u::AbstractVector,Δx::Float64,RHS;
        order=2,a=1.0,b=1.0,forcing=false)
    # Get penalties
    τ = SATpenalties(Robin,a,Δx,order)
    SAT = zeros(Float64,order)

    SAT[end] = τ * (b*u[end] - a*BDₓ(u,Right,Δx,order) + RHS[end])
    return SAT
end
