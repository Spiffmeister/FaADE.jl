

struct Injection_Dirichlet{
        TN<:NodeType,
        F1<:Function}
    side    :: TN
    RHS     :: F1
end



"""
    SAT_Dirichlet{
        TN<:NodeType,
        COORD,
        TT<:Real,
        VT<:Vector{TT},
        F1<:Function, LAT<:Function} <: SimultanousApproximationTerm{:Dirichlet}

Storage of all objects needed for a Dirichlet SAT ``u(x_i) - g(t) = 0``.

In Cartesian coordinates the SAT reads

``\\operatorname{SAT}_D = \\tau_0 H^{-1} E H^{-1} E (u - g) + \\tau_1 H^{-1} (K H D_x^T) H^{-1} E (u - g)``

In Curivlinear coordinates cross derivatives are included giving

``\\operatorname{SAT}_D = \\tau_0 H_q^{-1} E H_q^{-1} E (u - g) + \\tau_1 H^{-1}_q ((K_q H_q D_q^T) H_q^{-1} E + (K_{qr} H_r D_r^T) H_r^{-1} E) (u - g)``

where ``E`` picks out the relevant boundary term.

```
julia> RHS(X,t) = 0.0
julia> Δx = 1.0

julia> SAT_Dirichlet(RHS,Δx,Left,2)
```
"""
struct SAT_Dirichlet{
        TN<:NodeType,
        COORD,
        TT<:Real,
        VT<:Vector{TT},
        F1<:Function, LAT<:Function} <: SimultanousApproximationTerm{:Dirichlet}

    side        :: TN
    axis        :: Int
    order       :: Int
    RHS         :: F1
    H⁻¹EH⁻¹E    :: TT
    H⁻¹D₁ᵀE     :: VT
    Δx          :: TT
    τ₁          :: TT
    τ₀          :: VT
    loopaxis    :: LAT
    Δy          :: TT
    coordinates :: Symbol

    @doc """
        SAT_Dirichlet(RHS::F1,Δx::TT,side::NodeType{SIDE,AX},order::Int,Δy=TT(0),coord=:Cartesian;τ₀=nothing,τ₁=nothing) where {TT,F1,SIDE,AX}
    
    Constructor for the Dirichlet SAT.

    Required arguments:
    - `RHS`: Right hand side function
    - `Δx`: Grid spacing
    - `side`: Boundary side, either `Left`, `Right`, `Up`, `Down`
    - `order`: Order of the ``D_x`` operators, either `2` or `4`
    """   
    function SAT_Dirichlet(RHS::F1,Δx::TT,side::NodeType{SIDE,AX},order::Int,Δy=TT(0),coord=:Cartesian;τ₀=nothing,τ₁=nothing) where {TT,F1,SIDE,AX}
        # fullsat = "τH⁻¹ E H⁻¹E(u-f) + α H⁻¹ (K H D₁ᵀ) H⁻¹ E (u-f)"
    
        check_boundary(side)
    
        loopaxis = SelectLoopDirection(AX)
    
        if τ₀ === nothing
            τ₀ = TT(1)
        end
        if τ₁ === nothing
            τ₁ = [-TT(2)]
        end
    
        Hinv = _InverseMassMatrix(order,Δx,side)
        D₁ᵀ = _DerivativeTranspose(order,Δx,side)
        E = _BoundaryOperator(TT,side)
        
        # τ * H⁻¹ * E * SD.H⁻¹ * E * (u-f)
        if SIDE == :Left
            H⁻¹EH⁻¹E    = Hinv[1]*E*Hinv[1]*E
        elseif SIDE == :Right
            H⁻¹EH⁻¹E    = Hinv[end]*E*Hinv[end]*E
        end
        
        # α H⁻¹ * (K H D₁ᵀ) * H⁻¹ * E * (u-f)
        H⁻¹D₁ᵀE     = Hinv.*D₁ᵀ.*E
    
        new{typeof(side),coord,TT,Vector{TT},typeof(RHS),typeof(loopaxis)}(side,AX,order,RHS,H⁻¹EH⁻¹E,H⁻¹D₁ᵀE,Δx,τ₀,τ₁,loopaxis,Δy,coord)
    end
end





###
"""
    SAT_Dirichlet_explicit!

Dirichlet SAT for explicit solvers. Currently no explicit solvers are implemented so these haven't been tested and are not used anywhere.
"""
function SAT_Dirichlet_explicit! end
function SAT_Dirichlet_explicit!(dest::VT,u::VT,RHS::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:Union{NodeType{:Left},NodeType{:Up}}}
    for i = 1:SD.order
        dest[i] += SD.τ₁*c[1]*SD.H⁻¹D₁ᵀE[i]*(u[1] - RHS[1])
    end
    S[1] += SD.τ₀(c[1])*SD.H⁻¹EH⁻¹E*(u[1] - RHS[1])
end
function SAT_Dirichlet_explicit!(dest::VT,u::VT,RHS::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:Union{NodeType{:Right},NodeType{:Down}}}
    j = lastindex(dest)
    for i = 1:SD.order
        dest[j-SD.order+i] += SD.τ₁*C[j]*SD.H⁻¹D₁ᵀE[i]*(u[j] - RHS[j]) #D₁ᵀE₀
    end
    S[j] += SD.τ₀(c[j])*c[j]SD.H⁻¹EH⁻¹E*(u[j] - RHS[j])
end




# Operation along vector
"""
Dirichlet SAT method for implicit solvers. Applys boundary data i.e. the portion of the SAT which operates on ``g`` only.
    
All methods wrap around the 1D method.

    SAT_Dirichlet_data!(dest::VT,data::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:NodeType{SIDE}} where SIDE
        
The 2D method loops over the 1D method using the `loopaxis` value in [`SAT_Dirichlet`](@ref).
        
    SAT_Dirichlet_data!(dest::AT,data::AT,c::AT,SD::SAT_Dirichlet{TN,:Cartesian,TT}) where {AT<:AbstractMatrix,TN<:NodeType,TT}

The curvilinear call wraps around the 1D method and then calls the [`FirstDerivativeTranspose`](@ref) for the cross derivative term.

    SAT_Dirichlet_data!(dest::AT,data::AT,cx::KT,cxy::KT,SD::SAT_Dirichlet{TN,:Curvilinear,TT}) where {AT<:AbstractMatrix,TN<:NodeType,KT,TT}
"""
function SAT_Dirichlet_data! end
function SAT_Dirichlet_data!(dest::VT,data::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:NodeType{SIDE}} where SIDE

    # Pick out the correct index of the boundary element
    SIDE == :Left ? j = 1 : j = lastindex(dest) # First or last index
    SIDE == :Left ? m = 0 : m = j-SD.order      # Offset for looping
    SIDE == :Left ? o = 1 : o = lastindex(data) # Offset for data

    dest[j] += -SD.τ₀[1]*SD.H⁻¹EH⁻¹E*c[j]*data[o]
    for i = 1:SD.order #nodes SD.nodes
        dest[m+i] += -SD.τ₁*c[j]*SD.H⁻¹D₁ᵀE[i]*data[o] #u[Left]
    end
    dest
end
function SAT_Dirichlet_data!(dest::AT,data::AT,c::AT,SD::SAT_Dirichlet{TN,:Cartesian,TT}) where {AT<:AbstractMatrix,TN<:NodeType,TT}
    for (DEST,DATA,C) in zip(SD.loopaxis(dest),SD.loopaxis(data),SD.loopaxis(c))
        SAT_Dirichlet_data!(DEST,DATA,C,SD)
    end
    dest
end
function SAT_Dirichlet_data!(dest::AT,data::AT,cx::KT,cxy::KT,SD::SAT_Dirichlet{TN,:Curvilinear,TT}) where {AT<:AbstractMatrix,TN<:NodeType,KT,TT}
    for (DEST,DATA,C) in zip(SD.loopaxis(dest),SD.loopaxis(data),SD.loopaxis(cx))
        SAT_Dirichlet_data!(DEST,DATA,C,SD)
    end
    n = size(dest,SD.axis)
    m = size(dest,mod1(SD.axis+1,2))
    
    if SD.side == Left
        DEST =  view(dest,  1,1:m)
        SRC =   view(data,  1,1:m)
        C =     view(-SD.τ₁*cxy,  1,1:m)
    elseif SD.side == Right
        DEST =  view(dest,  n,1:m)
        SRC =   view(data,  1,1:m)
        C =     view(-SD.τ₁*cxy,  n,1:m)
    elseif SD.side == Down
        DEST =  view(dest,  1:m,1)
        SRC =   view(data,  1:m,1)
        C =     view(-SD.τ₁*cxy,  1:m,1)
    elseif SD.side == Up
        DEST =  view(dest,  1:m,n)
        SRC =   view(data,  1:m,1)
        C =     view(-SD.τ₁*cxy,  1:m,n)
    end

    FirstDerivativeTranspose!(DEST,SRC,C,m,SD.Δy,SD.order,TT(1))
    dest
end



"""
Dirichlet SAT method for implicit solvers. Applys portion of SAT related to the solution.

All methods wrap around the 1D method.

    SAT_Dirichlet_solution!(dest::VT,data::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:NodeType{SIDE}} where SIDE

The 2D method loops over the 1D method using the `loopaxis` value in [`SAT_Dirichlet`](@ref).

    SAT_Dirichlet_solution!(dest::AT,data::AT,c::KT,SD::SAT_Dirichlet{TN,:Cartesian,TT}) where {AT<:AbstractMatrix,TN<:NodeType,TT,KT}

The curvilinear call wraps around the 1D method and then calls the [`FirstDerivativeTranspose`](@ref) for the cross derivative term.

    SAT_Dirichlet_solution!(dest::AT,data::AT,cx::KT,cxy::KT,SD::SAT_Dirichlet{TN,:Curvilinear,TT}) where {AT<:AbstractMatrix,TN<:NodeType,TT,KT}
"""
function SAT_Dirichlet_solution! end
function SAT_Dirichlet_solution!(dest::VT,data::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:NodeType{SIDE}} where SIDE

    SIDE == :Left ? j = 1 : j = lastindex(dest)
    SIDE == :Left ? m = 0 : m = j-SD.order

    dest[j] += SD.τ₀[1]*SD.H⁻¹EH⁻¹E*c[j]*data[j]
    for i = 1:SD.order #nodes SD.nodes
        dest[m+i] += SD.τ₁*c[j]*SD.H⁻¹D₁ᵀE[i]*data[j] #u[Right]
    end
    dest
end
function SAT_Dirichlet_solution!(dest::AT,data::AT,c::KT,SD::SAT_Dirichlet{TN,:Cartesian,TT}) where {AT<:AbstractMatrix,TN<:NodeType,TT,KT}
    for (DEST,DATA,C) in zip(SD.loopaxis(dest),SD.loopaxis(data),SD.loopaxis(c))
        SAT_Dirichlet_solution!(DEST,DATA,C,SD)
    end
    dest
end
function SAT_Dirichlet_solution!(dest::AT,data::AT,cx::KT,cxy::KT,SD::SAT_Dirichlet{TN,:Curvilinear,TT}) where {AT<:AbstractMatrix,TN<:NodeType,TT,KT}
    for (DEST,DATA,C) in zip(SD.loopaxis(dest),SD.loopaxis(data),SD.loopaxis(cx))
        SAT_Dirichlet_solution!(DEST,DATA,C,SD)
    end

    n = size(dest,SD.axis)
    m = size(dest,mod1(SD.axis+1,2))
    
    if SD.side == Left
        DEST =  view(dest,      1,1:m)
        SRC =   view(data,      1,1:m)
        C =     view(SD.τ₁*cxy,  1,1:m)
    elseif SD.side == Right
        DEST =  view(dest,      n,1:m)
        SRC =   view(data,      n,1:m)
        C =     view(SD.τ₁*cxy,  n,1:m)
    elseif SD.side == Down
        DEST =  view(dest,      1:m,1)
        SRC =   view(data,      1:m,1)
        C =     view(SD.τ₁*cxy,  1:m,1)
    elseif SD.side == Up
        DEST =  view(dest,      1:m,n)
        SRC =   view(data,      1:m,n)
        C =     view(SD.τ₁*cxy,  1:m,n)
    end

    FirstDerivativeTranspose!(DEST,SRC,C,m,SD.Δy,SD.order,TT(1))
    dest
end


