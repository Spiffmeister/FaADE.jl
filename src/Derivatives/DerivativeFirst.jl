#=====================================#
#====== FIRST DERIVATIVE METHODS =====#
#=====================================#
# Author: Dean Muir, Kenneth Duru


"""
    FirstDerivativeInternal(u::AbstractVector{T},Δx::T,::NodeType{:Internal},order::Integer)
Single node 1D first derivative function.

See also [`FirstDerivativeInternal!`](@ref) for in place multi-node functions
"""
@inline function FirstDerivativeInternal(u::AbstractVector{T},Δx::T,::NodeType{:Internal},order::Integer) where T
    if order == 2
        j = 2
        return (u[j+1] - u[j-1])/(T(2)*Δx)
    elseif order == 4
        j = 3
        return (T(1/12)*u[i-2] - T(2/3)*u[i-1] + T(2/3)*u[i+1] - T(1/12)*u[i+2])/Δx
    elseif order == 6
        j = 4
        return (-T(1/60)*u[i-3] + T(3/20)*u[j-2] - T(3/4)*u[i-1] + T(3/4)*u[i+1] - T(3/20)*u[i+2] + T(1/60)*u[i+3])/Δx
    end
end




"""
    FirstDerivativeInternal
Single node 1D first derivative function.
"""
function FirstDerivativeInternal end
@inline function FirstDerivativeInternal(u::AT,Δx::T,::DerivativeOrder{2},i::Integer) where {T,AT<:AbstractVector{T}}
    # for i = 2:n-1
    @inbounds (u[i+1] - u[i-1])/(T(2)*Δx)
    # end
end
@inline function FirstDerivativeInternal(u::AT,Δx::T,::DerivativeOrder{4},i::Integer) where {T,AT<:AbstractVector{T}}
    # for i = 3:n-2
    @inbounds (T(1/12)*u[i-2] - T(2/3)*u[i-1] + T(2/3)*u[i+1] - T(1/12)*u[i+2])/Δx
    # end
end
@inline function FirstDerivativeInternal(u::AT,Δx::T,::DerivativeOrder{6},i::Integer) where {T,AT<:AbstractVector{T}}
    # for i = 4:n-3
    @inbounds (-T(1/60)*u[i-3] + T(3/20)*u[i-2] - T(3/4)*u[i-1] + T(3/4)*u[i+1] - T(3/20)*u[i+2] + T(1/60)*u[i+3])/Δx
    # end
end

"""
    FirstDerivativeInternal!(uₓ::AT,u::AT,Δx::T,n::Int,DO::DerivativeOrder{O},α::T)
In place first derivative function for internal nodes
"""
@inline function FirstDerivativeInternal!(uₓ::AT,u::AT,Δx::T,n::Int,DO::DerivativeOrder{O},α::T) where {T,AT<:AbstractVector{T},O}
    m = floor(Int,(O+2)/2)
    for i = m:n-m+1
        @inbounds uₓ[i] = α*uₓ[i] + FirstDerivativeInternal(u,Δx,DO,i)
    end
end
"""
    FirstDerivativeInternal!(uₓ::AT,c::AT,u::AT,Δx::T,n::Int,DO::DerivativeOrder,α::T)
In place first derivative function for internal nodes for computing ``DₓC Dₓu``
"""
@inline function FirstDerivativeInternal!(uₓ::AT,c::AT,u::AT,Δx::T,n::Int,DO::DerivativeOrder{O},α::T) where {T,AT<:AbstractVector{T},O}
    local ux :: T
    m = floor(Int,(O+2)/2)
    for i = m:n-m+1
        ux = FirstDerivativeInternal(u,Δx,DO,i)
        uₓ[i] = α*uₓ[i] + c[i]*ux
    end
end




"""
    FirstDerivativeBoundary!
1D in place function for first derivative on boundary nodes
"""
function FirstDerivativeBoundary! end
@inline function FirstDerivativeBoundary!(uₓ::AT,c::AT,u::AT,Δx::T,node::NodeType{TN},::DerivativeOrder{2},α::T) where {T,AT<:AbstractVector{T},TN}
    node == :Left ? i = 1 : i = -1
    node == :Left ? j = 1 : j = lastindex(u)
    uₓ[j] = α*uₓ[j] + c[j]*T(i)*(u[j+i] - u[j])/Δx
end
@inline function FirstDerivativeBoundary!(uₓ::AT,c::AT,u::AT,Δx::T,node::NodeType{TN},::DerivativeOrder{4},α::T) where {T,AT<:AbstractVector{T},TN}
    node == :Left ? i = 1 : i = -1
    node == :Left ? j = 1 : j = n
    uₓ[j]       = α*uₓ[i]   + c[j]      * T(i)*(T(24/17)*u[j] + T(59/34)*u[j+i]   - T(4/17)*u[j+2i]   - T(3/34)*u[j+3i])/Δx
    uₓ[j+i]     = α*uₓ[j+i] + c[j+i]    * T(i)*(T(1/2)*u[j]   + T(1/2)*u[j+2i])/Δx
    uₓ[j+2i]    = α*uₓ[j+2i]+ c[j+2i]   * T(i)*(T(4/43)*u[j]  - T(59/86)*u[j+i]   + T(59/86)*u[j+3i]  - T(4/43)*u[j+4i])/Δx
    uₓ[j+3i]    = α*uₓ[j+2i]+ c[j+3i]   * T(i)*(T(3/98)*u[j]  - T(59/98)*u[j+2i]  + T(32/49)*u[j+4i]  - T(4/49)*u[j+5i])/Δx
    # uₓ[j:i:j+3i] = uₓ[j:i:j+i]/Δx
end
@inline function FirstDerivativeBoundary!(uₓ::AT,c::AT,u::AT,Δx::T,node::TN,::DerivativeOrder{6},α::T) where {T,AT<:AbstractVector{T},TN<:Union{NodeType{:Left},NodeType{:Right}}}
    node <: Left ? i = 1 : i = -1
    node <: Left ? j = 1 : j = n
    uₓ[j]       = α*uₓ[j]   + c[j]      *T(i)*( T(-1.582533518939116)*u[j] + T(2.033378678700676)*u[j+i] - T(0.141512858744873)*u[j+2i] + T(-0.450398306578272)*u[j+3i] + T(0.104488069284042)*u[j+4i] + T(0.036577936277544)*u[j+5i] )/Δx
    uₓ[j+i]     = α*uₓ[j+i] + c[j+i]    *T(i)*( T(-0.462059195631158)*u[j] + T(0.287258622978251)*u[j+2i] + T(0.258816087376832)*u[j+3i] + T(-0.069112065532624)*u[j+4i] - T(0.014903449191300)*u[j+5i] )/Δx
    uₓ[j+2i]    = α*uₓ[j+2i]+ c[j+2i]   *T(i)*( T(0.071247104721830)* u[j] - T(0.636451095137907)*u[j+i] + T(0.606235523609147)*u[j+3i] + T(-0.022902190275815)*u[j+4i] - T(0.018129342917256)*u[j+5i] )/Δx
    uₓ[j+3i]    = α*uₓ[j+3i]+ c[j+3i]   *T(i)*( T(0.114713313798970)* u[j] - T(0.290087484386815)*u[j+i] - T(0.306681191361148)*u[j+2i] + T(0.520262285050482)*u[j+4i]  - T(0.051642265516119)*u[j+5i] + T(0.013435342414630)*u[j+6i] )/Δx
    uₓ[j+4i]    = α*uₓ[j+4i]+ c[j+4i]   *T(i)*( T(-0.036210680656541)*u[j] + T(0.105400944933782)*u[j+i] + T(0.015764336127392)*u[j+2i] + T(-0.707905442575989)*u[j+3i] + T(0.769199413962647)*u[j+5i] - T(0.164529643265203)*u[j+6i] + T(0.018281071473911)*u[j+7i] )/Δx
    uₓ[j+5i]    = α*uₓ[j+5i]+ c[j+5i]   *T(i)*( T(-0.011398193015050)*u[j] + T(0.020437334208704)*u[j+i] + T(0.011220896474665)*u[j+2i] + T( 0.063183694641876)*u[j+3i] - T(0.691649024426814)*u[j+4i] + T(0.739709139060752)*u[j+6i] + T(-0.147941827812150)*u[j+7i] + T(0.016437980868017)*u[j+8i] )/Δx
    # uₓ[j:i:j+5i] = uₓ[j:i:j+5i]/Δx
end

@inline function FirstDerivativeBoundary!(uₓ::AT,u::AT,Δx::T,node::NodeType{TN},::DerivativeOrder{2},α::T) where {T,AT<:AbstractVector{T},TN}
    node == :Left ? i = 1 : i = -1
    node == :Left ? j = 1 : j = lastindex(u)
    uₓ[j] = α*uₓ[j] + T(i)*(u[j+i] - u[j])/Δx
end
@inline function FirstDerivativeBoundary!(uₓ::AT,u::AT,Δx::T,node::TN,::DerivativeOrder{4},α::T) where {T,AT<:AbstractVector{T},TN<:Union{NodeType{:Left},NodeType{:Right}}}
    node == Left ? i = 1 : i = -1
    node == Left ? j = 1 : j = lastindex(u)
    uₓ[j]       = T(i)*(T(-24/17)*u[j]  + T(59/34)*u[j+i]       + T(-4/17)*u[j+2i] + T(-3/34)*u[j+3i])/Δx
    uₓ[j+i]     = T(i)*(T(-1/2)*u[j]    + T(1/2)*u[j+2i])/Δx
    uₓ[j+2i]    = T(i)*(T(4/43)*u[j]    + T(-59/86)*u[j+i]      + T(59/86)*u[j+3i] + T(-4/43)*u[j+4i])/Δx
    uₓ[j+3i]    = T(i)*(T(3/98)*u[j]    + T(-59/98)*u[j+2i]     + T(32/49)*u[j+4i] + T(-4/49)*u[j+5i])/Δx
end
@inline function FirstDerivativeBoundary!(uₓ::AT,u::AT,Δx::T,node::TN,::DerivativeOrder{6},α::T) where {T,AT<:AbstractVector{T},TN<:Union{NodeType{:Left},NodeType{:Right}}}
    node == Left ? i = 1 : i = -1
    node == Left ? j = 1 : j = n
    uₓ[j]       = α*uₓ[j]   + T(i)*( T(-1.582533518939116)*u[j] + T(2.033378678700676)*u[j+i] - T(0.141512858744873)*u[j+2i] + T(-0.450398306578272)*u[j+3i] + T(0.104488069284042)*u[j+4i] + T(0.036577936277544)*u[j+5i] )/Δx
    uₓ[j+i]     = α*uₓ[j+i] + T(i)*( T(-0.462059195631158)*u[j] + T(0.287258622978251)*u[j+2i] + T(0.258816087376832)*u[j+3i] + T(-0.069112065532624)*u[j+4i] - T(0.014903449191300)*u[j+5i] )/Δx
    uₓ[j+2i]    = α*uₓ[j+2i]+ T(i)*( T(0.071247104721830)* u[j] - T(0.636451095137907)*u[j+i] + T(0.606235523609147)*u[j+3i] + T(-0.022902190275815)*u[j+4i] - T(0.018129342917256)*u[j+5i] )/Δx
    uₓ[j+3i]    = α*uₓ[j+3i]+ T(i)*( T(0.114713313798970)* u[j] - T(0.290087484386815)*u[j+i] - T(0.306681191361148)*u[j+2i] + T(0.520262285050482)*u[j+4i]  - T(0.051642265516119)*u[j+5i] + T(0.013435342414630)*u[j+6i] )/Δx
    uₓ[j+4i]    = α*uₓ[j+4i]+ T(i)*( T(-0.036210680656541)*u[j] + T(0.105400944933782)*u[j+i] + T(0.015764336127392)*u[j+2i] + T(-0.707905442575989)*u[j+3i] + T(0.769199413962647)*u[j+5i] - T(0.164529643265203)*u[j+6i] + T(0.018281071473911)*u[j+7i] )/Δx
    uₓ[j+5i]    = α*uₓ[j+5i]+ T(i)*( T(-0.011398193015050)*u[j] + T(0.020437334208704)*u[j+i] + T(0.011220896474665)*u[j+2i] + T( 0.063183694641876)*u[j+3i] - T(0.691649024426814)*u[j+4i] + T(0.739709139060752)*u[j+6i] + T(-0.147941827812150)*u[j+7i] + T(0.016437980868017)*u[j+8i] )/Δx
    # uₓ[j:i:j+5i] = uₓ[j:i:j+5i]/Δx
end





