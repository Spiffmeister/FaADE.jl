


"""
    SAT_Neumann{TN<:NodeType,COORD,TT<:Real,VT<:Vector{TT},F1<:Function, LAT<:Function} <: SimultanousApproximationTerm{:Neumann}

Simulatenous approximation term for Neumann boundary conditions 
    ``\\left.\\frac{\\partial u}{\\partial x}\\right|_{x_i} = g(t) \\iff \\frac{\\partial}{\\partial x} u(x_i) - g(t) = 0``

In Cartesian coordinates the SAT reads

``\\operatorname{SAT}_N = \\tau H^{-1} E (K D_x u - g)``

In curvilinear coordinates the SAT is

``\\operatorname{SAT}_N = \\tau H_x^{-1} E (K_x D_x u - g) + \\tau H_x^{-1} E K_{xy} D_y u``
"""
struct SAT_Neumann{
    TN<:NodeType,
    COORD,
    TT<:Real,
    VT<:Vector{TT},
    F1<:Function, LAT<:Function} <: SimultanousApproximationTerm{:Neumann}
    
    side        :: TN
    axis        :: Int
    order       :: Int
    RHS         :: F1
    H⁻¹E        :: TT
    D₁          :: VT
    Δx          :: TT
    τ           :: TT
    loopaxis    :: LAT
    Δy          :: TT
    coordinates :: Symbol
    ### CONSTRUCTOR ###

    @doc """
        SAT_Neumann(RHS::F1,Δx::TT,side::NodeType{SIDE,AX},order::Int,Δy=TT(0),coord=:Cartesian;τ=nothing,Δy=TT(0),coord=:Cartesian) where {TT,F1,SIDE,AX}
    
    Constructor for the Neumann SAT

    Required arguments:
    - `RHS`: Right hand side function
    - `Δx`: Grid spacing
    - `side`: Boundary side, either `Left`, `Right`, `Up`, `Down`
    - `order`: Order of the ``D_x`` operators, either `2` or `4`
    """
    function SAT_Neumann(RHS::F1,Δx::TT,side::NodeType{SIDE,AX},order::Int,Δy=TT(0),coord=:Cartesian;τ=nothing) where {TT,F1,SIDE,AX}

        check_boundary(side)

        LA = SelectLoopDirection(AX)


        Hinv = _InverseMassMatrix(order,Δx,side)
        # E = _BoundaryOperator(TT,side)
        D₁ = _BoundaryDerivative(order,Δx,side)


        # τ * H⁻¹ * E
        if SIDE == :Left
            H⁻¹E    = Hinv[1]
        elseif SIDE == :Right
            n = -1 # normal vector
            H⁻¹E    = n*Hinv[end]
        end

        if isnothing(τ)
            τ = TT(1)
        end

        new{typeof(side),coord,TT,Vector{TT},F1,typeof(LA)}(side,AX,order,RHS,H⁻¹E,D₁,Δx,τ,LA,Δy,coord)
    end
end




#=== Explicit methods ===#
"""
    SAT_Neumann_explicit!
Neumann boundary SAT for explicit solvers. No explicit solvers are implemented so this has not been tested and is not used.
"""
function SAT_Neumann_explicit! end
function SAT_Neumann_explicit!(RHS::Function,SAT::AbstractArray{T},::NodeType{:Left},u::AbstractArray{T},c::AbstractArray{T},t::Float64,
        τ::T,BD::AbstractArray{T},
        order::Int,loopaxis::Function) where T

    for (S,C,U) in zip(loopaxis(SAT),loopaxis(u),loopaxis(c))
        S[1] = τ*(C[1] * dot(BD,U[1:order]) - RHS(t))
    end
end
function SAT_Neumann_explicit!(RHS::Function,SAT::AbstractArray{T},::NodeType{:Right},u::AbstractArray{T},c::AbstractArray{T},t::Float64,
        τ::T,BD::AbstractArray{T},
        order::Int,loopaxis::Function) where T

    for (S,C,U) in zip(loopaxis(SAT),loopaxis(u),loopaxis(c))
        S[end] -= τ*(C[end] * dot(BD,U[end-order+1:end]) - RHS(t))
    end
end







#== NEW ==#
"""
Neumann SAT method for implicit solvers. Applies boundary data i.e. the portion of the SAT which operates on ``g`` only.

All methods wrap around the 1D method.

    SAT_Neumann_data!(dest::VT,u::VT,SN::SAT_Neumann{TN}) where {VT<:AbstractVector,TN<:NodeType{SIDE}} where SIDE

The 2D method loops over the 1D method using the `loopaxis` value in [`SAT_Neumann`](@ref).

    SAT_Neumann_data!(dest::AT,u,SN::SAT_Neumann{TN,:Cartesian,TT}) where {AT<:AbstractMatrix,TN,TT}

The curvilinear call wraps around the 1D method.

    SAT_Neumann_data!(dest::AT,data::AT,SN::SAT_Neumann{TN,:Curvilinear,TT}) where {AT<:AbstractMatrix,TN<:NodeType,TT}

"""
function SAT_Neumann_data! end
function SAT_Neumann_data!(dest::VT,u::VT,SN::SAT_Neumann{TN}) where {VT<:AbstractVector,TN<:NodeType{SIDE}} where SIDE

    SIDE == :Left ? j = 1 : j = lastindex(dest)
    SIDE == :Left ? o = 1 : o = lastindex(u)
    
    dest[j] -= SN.τ*SN.H⁻¹E*u[o]
end
function SAT_Neumann_data!(dest::AT,u,SN::SAT_Neumann{TN,:Cartesian,TT}) where {AT<:AbstractMatrix,TN,TT}
    for (DEST,U) in zip(SN.loopaxis(dest),SN.loopaxis(u))
        SAT_Neumann_data!(DEST,U,SN)
    end
end
function SAT_Neumann_data!(dest::AT,data::AT,SN::SAT_Neumann{TN,:Curvilinear,TT}) where {AT<:AbstractMatrix,TN<:NodeType,TT}
    for (DEST,DATA) in zip(SN.loopaxis(dest),SN.loopaxis(data))
        SAT_Neumann_data!(DEST,DATA,SN)
    end
    dest
end



"""
Neumann SAT method for implicit solvers. Applys portion of SAT related to the solution.

All methods wrap around the 1D method.

    SAT_Neumann_solution!(dest::VT,u::VT,c::VT,SN::SAT_Neumann{TN}) where {VT<:AbstractVector,TN<:NodeType{SIDE}} where SIDE

The 2D method loops over the 1D method using the `loopaxis` value in [`SAT_Neumann`](@ref).

    SAT_Neumann_solution!(dest::AT,u::AT,c::AT,SN::SAT_Neumann{TN,:Cartesian,TT}) where {AT<:AbstractMatrix,TN,TT}

The curvilinear call wraps around the 1D method. Calls the [`D₁!`](@ref) operator from [`FaADE.Derivatives`](@ref)

    SAT_Neumann_solution!(dest::AT,data::AT,c::AT,cxy::AT,SN::SAT_Neumann{TN,:Curvilinear,TT}) where {AT,TN<:NodeType,TT}

"""
function SAT_Neumann_solution! end
function SAT_Neumann_solution!(dest::VT,u::VT,c::VT,SN::SAT_Neumann{TN}) where {VT<:AbstractVector,TN<:NodeType{SIDE}} where SIDE

    SIDE == :Left ? j = 1 : j = lastindex(dest)
    SIDE == :Left ? m = 0 : m = j - SN.order

    for i = 1:SN.order
        dest[j] += SN.τ * SN.H⁻¹E * c[j] * SN.D₁[i]*u[m+i]
    end
end
function SAT_Neumann_solution!(dest::AT,u::AT,c::AT,SN::SAT_Neumann{TN,:Cartesian,TT}) where {AT<:AbstractMatrix,TN,TT}
    for (DEST,U,C) in zip(SN.loopaxis(dest),SN.loopaxis(u),SN.loopaxis(c))
        SAT_Neumann_solution!(DEST,U,C,SN)
    end
end
function SAT_Neumann_solution!(dest::AT,data::AT,c::AT,cxy::AT,SN::SAT_Neumann{TN,:Curvilinear,TT}) where {AT,TN<:NodeType,TT}
    # D_q u
    for (DEST,DATA,C) in zip(SN.loopaxis(dest),SN.loopaxis(data),SN.loopaxis(c))
        SAT_Neumann_solution!(DEST,DATA,C,SN)
    end
    n = size(dest,SN.axis)
    m = size(dest,mod1(SN.axis+1,2))
    
    #Compute D_r u
    if SN.side == Left
        DEST =  view(dest,              1,1:m)
        SRC =   view(data,              1,1:m)
        C =     view(SN.τ*SN.H⁻¹E*cxy,  1,1:m)
    elseif SN.side == Right
        DEST =  view(dest,              n,1:m)
        SRC =   view(data,              1,1:m)
        C =     view(SN.τ*SN.H⁻¹E*cxy,  n,1:m)
    elseif SN.side == Down
        DEST =  view(dest,              1:m,1)
        SRC =   view(data,              1:m,1)
        C =     view(SN.τ*SN.H⁻¹E*cxy,  1:m,1)
    elseif SN.side == Up
        DEST =  view(dest,              1:m,n)
        SRC =   view(data,              1:m,1)
        C =     view(SN.τ*SN.H⁻¹E*cxy,  1:m,n)
    end
    D₁!(DEST,C,SRC,m,SN.Δy,Val(SN.order),TT(1))
    dest
end

