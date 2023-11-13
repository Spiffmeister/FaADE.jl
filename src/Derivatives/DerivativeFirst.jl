#=====================================#
#====== FIRST DERIVATIVE METHODS =====#
#=====================================#
# Author: Dean Muir, Kenneth Duru


"""
    FirstDerivativeInternal
Single node 1D first derivative function.
"""
function FirstDerivativeInternal end
@inline function FirstDerivativeInternal(u::AT,Δx::T,::DerivativeOrder{2},i::Integer,β::T) where {T,AT<:AbstractVector{T}}
    @inbounds β*(u[i+1] - u[i-1])/(T(2)*Δx)
end
@inline function FirstDerivativeInternal(u::AT,Δx::T,::DerivativeOrder{4},i::Integer,β::T) where {T,AT<:AbstractVector{T}}
    @inbounds β*(T(1/12)*u[i-2] - T(2/3)*u[i-1] + T(2/3)*u[i+1] - T(1/12)*u[i+2])/Δx
end
@inline function FirstDerivativeInternal(u::AT,Δx::T,::DerivativeOrder{6},i::Integer,β::T) where {T,AT<:AbstractVector{T}}
    @inbounds β*(-T(1/60)*u[i-3] + T(3/20)*u[i-2] - T(3/4)*u[i-1] + T(3/4)*u[i+1] - T(3/20)*u[i+2] + T(1/60)*u[i+3])/Δx
end


"""
    FirstDerivativeInternal!(uₓ::AT,u::AT,Δx::T,n::Int,DO::DerivativeOrder{O},α::T)
In place first derivative function for internal nodes
"""
function FirstDerivativeInternal! end
@inline function FirstDerivativeInternal!(uₓ::AT,u::AT,Δx::TT,n::Int,DO::DerivativeOrder{O},α::TT) where {TT,AT<:AbstractVector{TT},O}
    m = floor(Int,(O+2)/2)
    for i = m:n-m+1
        @inbounds uₓ[i] = α*uₓ[i] + FirstDerivativeInternal(u,Δx,DO,i,TT(1))
    end
end
@inline function FirstDerivativeInternal!(uₓ::AT,K::AT,u::AT,Δx::TT,n::Int,DO::DerivativeOrder{O},α::TT) where {TT,AT<:AbstractVector{TT},O}
    m = floor(Int,(O+2)/2)
    for i = m:n-m+1
        @inbounds uₓ[i] = α*uₓ[i] + FirstDerivativeInternal(u,Δx,DO,i,K[i]) 
    end
end


"""
    FirstDerivativeBoundary!
1D in place function for first derivative on boundary nodes
"""
function FirstDerivativeBoundary! end
@inline function FirstDerivativeBoundary!(uₓ::AT,u::AT,Δx::TT,NT::NodeType,::DerivativeOrder{2},α::TT) where {TT,AT<:AbstractVector{TT}}
    NT == Left ? i = 1 : i = -1
    NT == Left ? j = 1 : j = length(u)
    uₓ[j] = α*uₓ[j] + (u[j+i] - u[j])/Δx
end
@inline function FirstDerivativeBoundary!(uₓ::AT,u::AT,Δx::TT,NT::NodeType,::DerivativeOrder{4},α::TT) where {TT,AT<:AbstractVector{TT}}
    NT == Left ? i = 1 : i = -1
    NT == Left ? j = 1 : j = length(u)
    uₓ[j]       = α*uₓ[j]   + TT(i)*(TT(-24/17)*u[j] + TT(59/34)*u[j+i]   + TT(-4/17)*u[j+2i]   + TT(-3/34)*u[j+3i])/Δx
    uₓ[j+i]     = α*uₓ[j+i] + TT(i)*(TT(-1/2)*u[j]   + TT(1/2)*u[j+2i])/Δx
    uₓ[j+2i]    = α*uₓ[j+2i]+ TT(i)*(TT(4/43)*u[j]  + TT(-59/86)*u[j+i]   + TT(59/86)*u[j+3i]  + TT(-4/43)*u[j+4i])/Δx
    uₓ[j+3i]    = α*uₓ[j+3i]+ TT(i)*(TT(3/98)*u[j]  + TT(-59/98)*u[j+2i]  + TT(32/49)*u[j+4i]  + TT(-4/49)*u[j+5i])/Δx
end
@inline function FirstDerivativeBoundary!(uₓ::AT,u::AT,Δx::TT,NT::NodeType,::DerivativeOrder{6},α::TT) where {TT,AT<:AbstractVector{TT}}
    NT == Left ? i = 1 : i = -1
    NT == Left ? j = 1 : j = length(u)
    uₓ[j]       = α*uₓ[j]   + TT(i)*( TT(-1.582533518939116)*u[j] + TT(2.033378678700676)*u[j+i] - TT(0.141512858744873)*u[j+2i] + TT(-0.450398306578272)*u[j+3i] + TT(0.104488069284042)*u[j+4i] + TT(0.036577936277544)*u[j+5i] )/Δx
    uₓ[j+i]     = α*uₓ[j+i] + TT(i)*( TT(-0.462059195631158)*u[j] + TT(0.287258622978251)*u[j+2i] +TT(0.258816087376832)*u[j+3i] + TT(-0.069112065532624)*u[j+4i] - TT(0.014903449191300)*u[j+5i] )/Δx
    uₓ[j+2i]    = α*uₓ[j+2i]+ TT(i)*( TT(0.071247104721830)* u[j] - TT(0.636451095137907)*u[j+i] + TT(0.606235523609147)*u[j+3i] + TT(-0.022902190275815)*u[j+4i] - TT(0.018129342917256)*u[j+5i] )/Δx
    uₓ[j+3i]    = α*uₓ[j+3i]+ TT(i)*( TT(0.114713313798970)* u[j] - TT(0.290087484386815)*u[j+i] - TT(0.306681191361148)*u[j+2i] + TT(0.520262285050482)*u[j+4i]  - TT(0.051642265516119)*u[j+5i] + TT(0.013435342414630)*u[j+6i] )/Δx
    uₓ[j+4i]    = α*uₓ[j+4i]+ TT(i)*( TT(-0.036210680656541)*u[j] + TT(0.105400944933782)*u[j+i] + TT(0.015764336127392)*u[j+2i] + TT(-0.707905442575989)*u[j+3i] + TT(0.769199413962647)*u[j+5i] - TT(0.164529643265203)*u[j+6i] + T(0.018281071473911)*u[j+7i] )/Δx
    uₓ[j+5i]    = α*uₓ[j+5i]+ TT(i)*( TT(-0.011398193015050)*u[j] + TT(0.020437334208704)*u[j+i] + TT(0.011220896474665)*u[j+2i] + TT( 0.063183694641876)*u[j+3i] - TT(0.691649024426814)*u[j+4i] + TT(0.739709139060752)*u[j+6i] + TT(-0.147941827812150)*u[j+7i] + TT(0.016437980868017)*u[j+8i] )/Δx
end




@inline function FirstDerivativeBoundary!(uₓ::AT,K::AT,u::AT,Δx::TT,NT::NodeType,::DerivativeOrder{2},α::TT) where {TT,AT<:AbstractVector{TT}}
    NT == Left ? i = 1 : i = -1
    NT == Left ? j = 1 : j = length(u)
    uₓ[j] = α*uₓ[j] + K[i]*(u[j+i] - u[j])/Δx
end
@inline function FirstDerivativeBoundary!(uₓ::AT,K::AT,u::AT,Δx::TT,NT::NodeType,::DerivativeOrder{4},α::TT) where {TT,AT<:AbstractVector{TT}}
    NT == Left ? i = 1 : i = -1
    NT == Left ? j = 1 : j = length(u)
    uₓ[j]       = α*uₓ[j]   + K[j]      *TT(i)*(TT(-24/17)*u[j] + TT(59/34)*u[j+i]   + TT(-4/17)*u[j+2i]   + TT(-3/34)*u[j+3i])/Δx
    uₓ[j+i]     = α*uₓ[j+i] + K[j+i]    *TT(i)*(TT(-1/2)*u[j]   + TT(1/2)*u[j+2i])/Δx
    uₓ[j+2i]    = α*uₓ[j+2i]+ K[j+2i]   *TT(i)*(TT(4/43)*u[j]  + TT(-59/86)*u[j+i]   + TT(59/86)*u[j+3i]  + TT(-4/43)*u[j+4i])/Δx
    uₓ[j+3i]    = α*uₓ[j+3i]+ K[j+3i]   *TT(i)*(TT(3/98)*u[j]  + TT(-59/98)*u[j+2i]  + TT(32/49)*u[j+4i]  + TT(-4/49)*u[j+5i])/Δx
end
@inline function FirstDerivativeBoundary!(uₓ::AT,K::AT,u::AT,Δx::TT,NT::NodeType,::DerivativeOrder{6},α::TT) where {TT,AT<:AbstractVector{TT}}
    NT == Left ? i = 1 : i = -1
    NT == Left ? j = 1 : j = length(u)
    uₓ[j]       = α*uₓ[j]   + K[j]      *TT(i)*( TT(-1.582533518939116)*u[j] + TT(2.033378678700676)*u[j+i] - TT(0.141512858744873)*u[j+2i] + TT(-0.450398306578272)*u[j+3i] + TT(0.104488069284042)*u[j+4i] + TT(0.036577936277544)*u[j+5i] )/Δx
    uₓ[j+i]     = α*uₓ[j+i] + K[j+i]    *TT(i)*( TT(-0.462059195631158)*u[j] + TT(0.287258622978251)*u[j+2i] +TT(0.258816087376832)*u[j+3i] + TT(-0.069112065532624)*u[j+4i] - TT(0.014903449191300)*u[j+5i] )/Δx
    uₓ[j+2i]    = α*uₓ[j+2i]+ K[j+2i]   *TT(i)*( TT(0.071247104721830)* u[j] - TT(0.636451095137907)*u[j+i] + TT(0.606235523609147)*u[j+3i] + TT(-0.022902190275815)*u[j+4i] - TT(0.018129342917256)*u[j+5i] )/Δx
    uₓ[j+3i]    = α*uₓ[j+3i]+ K[j+3i]   *TT(i)*( TT(0.114713313798970)* u[j] - TT(0.290087484386815)*u[j+i] - TT(0.306681191361148)*u[j+2i] + TT(0.520262285050482)*u[j+4i]  - TT(0.051642265516119)*u[j+5i] + TT(0.013435342414630)*u[j+6i] )/Δx
    uₓ[j+4i]    = α*uₓ[j+4i]+ K[j+4i]   *TT(i)*( TT(-0.036210680656541)*u[j] + TT(0.105400944933782)*u[j+i] + TT(0.015764336127392)*u[j+2i] + TT(-0.707905442575989)*u[j+3i] + TT(0.769199413962647)*u[j+5i] - TT(0.164529643265203)*u[j+6i] + T(0.018281071473911)*u[j+7i] )/Δx
    uₓ[j+5i]    = α*uₓ[j+5i]+ K[j+5i]   *TT(i)*( TT(-0.011398193015050)*u[j] + TT(0.020437334208704)*u[j+i] + TT(0.011220896474665)*u[j+2i] + TT( 0.063183694641876)*u[j+3i] - TT(0.691649024426814)*u[j+4i] + TT(0.739709139060752)*u[j+6i] + TT(-0.147941827812150)*u[j+7i] + TT(0.016437980868017)*u[j+8i] )/Δx
end





@inline function FirstDerivativePeriodic(u::AT,Δx::T,::DerivativeOrder{2},n::Integer,i::Integer,β::T) where {T,AT<:AbstractVector{T}}
    @inbounds β*(u[_prev(i,n)] - u[_next(i,n)])/(T(2)*Δx)
end
@inline function FirstDerivativePeriodic(u::AT,Δx::T,::DerivativeOrder{4},n::Integer,i::Integer,β::T) where {T,AT<:AbstractVector{T}}
    @inbounds β*(T(1/12)*u[_prev(i-1,n)] - T(2/3)*u[_prev(i,n)] + T(2/3)*u[_next(i,n)] - T(1/12)*u[_next(i+1,n)])/Δx
end
@inline function FirstDerivativePeriodic(u::AT,Δx::T,::DerivativeOrder{6},n::Integer,i::Integer,β::T) where {T,AT<:AbstractVector{T}}
    @inbounds β*(-T(1/60)*u[_prev(i-2,n)] + T(3/20)*u[_prev(i-1,n)] - T(3/4)*u[_prev(i,n)] + T(3/4)*u[_next(i,n)] - T(3/20)*u[_next(i+1,n)] + T(1/60)*u[_next(i+2,n)])/Δx
end
function FirstDerivativePeriodic!(dest::VT,K::VT,u::VT,Δx::TT,DO::DerivativeOrder,n::Int,α::TT) where {TT,VT<:AbstractVector{TT}}
    for i = 1:n
        @inbounds dest[i] = α*dest[i] + FirstDerivativePeriodic(u,Δx,DO,n,i,K[i])
    end
end