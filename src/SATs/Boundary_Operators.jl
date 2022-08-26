"""
    BDₓᵀ

1. BDₓᵀ(u::AbstractVector,::NodeType{:Left},Δx::Float64,order::Int64=2)
2. BDₓᵀ(u::AbstractVector,::NodeType{:Right},Δx::Float64,order::Int64=2)

Transpose of first order derivative operator at the boundary needed to compute certain SATs
"""
function BDₓᵀ end
@inline function BDₓᵀ(u::Union{AbstractFloat,AbstractVector},::NodeType{:Left},Δx::Float64,order::Int64=2)
    if order == 2
        BOp = [-1.0, 1.0]
        BOp .*= u[1]/Δx
        return BOp
    elseif order == 4
        BOp = [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]
        BOp .*= u[1]/Δx
        return BOp
    elseif order == 6
        BOp = [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]
        BOp .*= u[1]/Δx
        return BOp
    else
        error("Order must be 2, 4, or 6")
    end
end
@inline function BDₓᵀ(u::Union{AbstractFloat,AbstractVector},::NodeType{:Right},Δx::Float64,order::Int64=2)
    if order == 2
        BOp = [-1.0, 1.0]
        BOp .*= u[end]/Δx
        return BOp
    elseif order == 4
        BOp = [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]
        BOp .*= u[end]/Δx
        return BOp
    elseif order == 6
        BOp = [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]
        BOp .*= u[end]/Δx
        return BOp
    else
        error("Order must be 2, 4, or 6")
    end
end





"""
    BDₓ

1. BDₓ(u::AbstractVector,::NodeType{:Left},Δx::Float64,order::Int64=2)
2. BDₓ(u::AbstractVector,::NodeType{:Right},Δx::Float64,order::Int64=2)

Transpose of first order derivative operator at the boundary needed to compute certain SATs
"""
function BDₓ end
@inline function BDₓ(u::AbstractVector,::NodeType{:Left},Δx::Float64,order::Int64=2)
    if order == 2
        #Boundary for the second order case
        return (u[2] - u[1])/Δx
    elseif order == 4
        #Boundary for the fourth order case
        B0p = [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]
        return sum(B0p.*u[1:4])/Δx
    elseif order == 6
        #Boundary for the sixth order case
        B0p = [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]
        return sum(B0p.*u[1:6])/Δx
    end
    return uₓ
end
@inline function BDₓ(u::AbstractVector,::NodeType{:Right},Δx::Float64,order::Int64=2)
    if order == 2
        #Boundary for the second order case
        return (u[end] - u[end-1])/Δx
    elseif order == 4
        #Boundary for the fourth order case
        Bnp = [3.0/34.0, 4.0/17.0, -59.0/34.0, 24.0/17.0]
        return sum(Bnp.*u[end-3:end])/Δx
    elseif order == 6
        #Boundary for the sixth order case
        Bnp = [-0.036577936277544, -0.104488069284042, 0.450398306578272, 0.141512858744873, -2.033378678700676, 1.582533518939116]
        return sum(Bnp.*u[end-5:end])/Δx
    end
end




"""
    SATpenalties

1. `SATpenalties(::BoundaryCondition{:Dirichlet,:Neumann},Δx::Float64,order::Int64)`
2. `SATpenalties(::BoundaryCondition{:Robin},a,Δx::Float64,order::Int64)`

Determines the penatly parameters for the given boundary conditions.
"""
function SATpenalties end
@inline function SATpenalties(::BoundaryCondition{:Dirichlet},Δx::Float64,order::Int64)
    # For reading in penalty parameters for Dirichlet SATs
    h = hval(order)

    α = 1.0 * (h * Δx)^-1 # α*H^{-1}

    τ = 1.0
    τ = -(1.0 + τ) * (h * Δx)^-2 # τ*H^{-1}H^{-1}
    return α, τ
end
@inline function SATpenalties(::BoundaryCondition{:Neumann},Δx::Float64,order::Int64)
    # For reading in penalty parameters for Neumann SATs
    h = hval(order)

    τ = 1.0/(h * Δx) # τ*H^{-1}
    return τ
end
@inline function SATpenalties(::BoundaryCondition{:Robin},a,Δx::Float64,order::Int64)
    h = hval(order)

    τ = 1.0/(a * h * Δx) # τ=1/a H^{-1}
    return τ
end


"""
    hval(order::Int64)

Returns the value of ``h^{-1}`` for the penalties
"""
@inline function hval(order::Int64)
    if order == 2
        return 0.5
    elseif order == 4
        return 17.0/48.0
    elseif order == 6
        return 13649.0/43200.0
    end
end




