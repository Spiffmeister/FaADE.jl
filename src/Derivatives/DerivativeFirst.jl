#=====================================#
#====== FIRST DERIVATIVE METHODS =====#
#=====================================#
# Author: Dean Muir, Kenneth Duru



"""
    FirstDerivativeInternal
Single node 1D first derivative function.
"""
function FirstDerivativeInternal end
@inline function FirstDerivativeInternal(u::AT,Δx::T,::Val{2},i::Integer,β::T) where {T,AT<:AbstractVector{T}}
    @inbounds β*(u[i+1] - u[i-1])/(T(2)*Δx)
end
@inline function FirstDerivativeInternal(u::AT,Δx::T,::Val{4},i::Integer,β::T) where {T,AT<:AbstractVector{T}}
    @inbounds β*(T(1/12)*u[i-2] - T(2/3)*u[i-1] + T(2/3)*u[i+1] - T(1/12)*u[i+2])/Δx
end
@inline function FirstDerivativeInternal(u::AT,Δx::T,::Val{6},i::Integer,β::T) where {T,AT<:AbstractVector{T}}
    @inbounds β*(-T(1/60)*u[i-3] + T(3/20)*u[i-2] - T(3/4)*u[i-1] + T(3/4)*u[i+1] - T(3/20)*u[i+2] + T(1/60)*u[i+3])/Δx
end



"""
    FirstDerivativeInternal!(uₓ::AT,u::AT,Δx::T,n::Int,DO::Val{O},α::T)
In place first derivative function for internal nodes
"""
function FirstDerivativeInternal! end
@inline function FirstDerivativeInternal!(dest::VT,u::AT,Δx::TT,n::Int,DO::Val{O},α::TT) where {TT,AT<:AbstractVector{TT},VT<:AbstractVector{TT},O}
    O == 2 ? m = O : m = O+1
    for i = m:n-m+1
        @inbounds dest[i] = α*dest[i] + FirstDerivativeInternal(u,Δx,DO,i,TT(1))
    end
    dest
end
@inline function FirstDerivativeInternal!(dest::VT,K::AT,u::AT,Δx::TT,n::Int,DO::Val{O},α::TT) where {TT,AT<:AbstractVector{TT},VT<:AbstractVector{TT},O}
    O == 2 ? m = O : m = O+1
    # m = floor(Int,(O+2)/2)
    for i = m:n-m+1
        @inbounds dest[i] = α*dest[i] + FirstDerivativeInternal(u,Δx,DO,i,K[i])
    end
    dest
end



"""
    FirstDerivativePeriodic
Single node 1D first derivative periodic stencil
"""
@inline function FirstDerivativePeriodic(u::AT,Δx::T,::Val{2},n::Integer,i::Integer,β::T) where {T,AT<:AbstractVector{T}}
    @inbounds β*(u[_prev(i,n)] - u[_next(i,n)])/(T(2)*Δx)
end
@inline function FirstDerivativePeriodic(u::AT,Δx::T,::Val{4},n::Integer,i::Integer,β::T) where {T,AT<:AbstractVector{T}}
    @inbounds β*(T(1/12)*u[_prev(i-1,n)] - T(2/3)*u[_prev(i,n)] + T(2/3)*u[_next(i,n)] - T(1/12)*u[_next(i+1,n)])/Δx
end
@inline function FirstDerivativePeriodic(u::AT,Δx::T,::Val{6},n::Integer,i::Integer,β::T) where {T,AT<:AbstractVector{T}}
    @inbounds β*(-T(1/60)*u[_prev(i-2,n)] + T(3/20)*u[_prev(i-1,n)] - T(3/4)*u[_prev(i,n)] + T(3/4)*u[_next(i,n)] - T(3/20)*u[_next(i+1,n)] + T(1/60)*u[_next(i+2,n)])/Δx
end

@inline function FirstDerivativePeriodic!(dest::VT,K::VT,u::VT,Δx::TT,DO::Val,n::Int,α::TT) where {TT,VT<:AbstractVector{TT}}

    # for i = 2:n-1
    #     dest[i] = α*dest[i] + FirstDerivativeInternal(u,Δx,DO,i,K[i])
    # end

    # dest[1] = α*dest[1] + K[1]*(u[2] - u[n-1])/(TT(2)*Δx)
    # dest[n] = α*dest[n] + K[n]*(u[2] - u[n-1])/(TT(2)*Δx)

    for i = 1:n
        @inbounds dest[i] = α*dest[i] + FirstDerivativePeriodic(u,Δx,DO,n,i,K[i])
    end
    dest
end


"""
    FirstDerivativeBoundary!
1D in place function for first derivative on boundary nodes
"""
function FirstDerivativeBoundary! end
@inline function FirstDerivativeBoundary!(uₓ::AT,u::AT,Δx::TT,NT::NodeType,::Val{2},α::TT) where {TT,AT<:AbstractVector{TT}}
    NT == Left ? i = 1 : i = -1
    NT == Left ? j = 1 : j = length(u)
    uₓ[j] = α*uₓ[j] + TT(i)*(u[j+i] - u[j])/Δx
end
@inline function FirstDerivativeBoundary!(uₓ::AT,u::AT,Δx::TT,NT::NodeType,::Val{4},α::TT) where {TT,AT<:AbstractVector{TT}}
    NT == Left ? i = 1 : i = -1
    NT == Left ? j = 1 : j = length(u)
    uₓ[j]       = α*uₓ[j]   + TT(i)*(TT(-24/17)*u[j] + TT(59/34)*u[j+i]   + TT(-4/17)*u[j+2i]   + TT(-3/34)*u[j+3i])/Δx
    uₓ[j+i]     = α*uₓ[j+i] + TT(i)*(TT(-1/2)*u[j]   + TT(1/2)*u[j+2i])/Δx
    uₓ[j+2i]    = α*uₓ[j+2i]+ TT(i)*(TT(4/43)*u[j]  + TT(-59/86)*u[j+i]   + TT(59/86)*u[j+3i]  + TT(-4/43)*u[j+4i])/Δx
    uₓ[j+3i]    = α*uₓ[j+3i]+ TT(i)*(TT(3/98)*u[j]  + TT(-59/98)*u[j+2i]  + TT(32/49)*u[j+4i]  + TT(-4/49)*u[j+5i])/Δx
end
@inline function FirstDerivativeBoundary!(uₓ::AT,u::AT,Δx::TT,NT::NodeType,::Val{6},α::TT) where {TT,AT<:AbstractVector{TT}}
    NT == Left ? i = 1 : i = -1
    NT == Left ? j = 1 : j = length(u)
    uₓ[j]       = α*uₓ[j]   + TT(i)*( TT(-1.582533518939116)*u[j] + TT(2.033378678700676)*u[j+i] - TT(0.141512858744873)*u[j+2i] + TT(-0.450398306578272)*u[j+3i] + TT(0.104488069284042)*u[j+4i] + TT(0.036577936277544)*u[j+5i] )/Δx
    uₓ[j+i]     = α*uₓ[j+i] + TT(i)*( TT(-0.462059195631158)*u[j] + TT(0.287258622978251)*u[j+2i] +TT(0.258816087376832)*u[j+3i] + TT(-0.069112065532624)*u[j+4i] - TT(0.014903449191300)*u[j+5i] )/Δx
    uₓ[j+2i]    = α*uₓ[j+2i]+ TT(i)*( TT(0.071247104721830)* u[j] - TT(0.636451095137907)*u[j+i] + TT(0.606235523609147)*u[j+3i] + TT(-0.022902190275815)*u[j+4i] - TT(0.018129342917256)*u[j+5i] )/Δx
    uₓ[j+3i]    = α*uₓ[j+3i]+ TT(i)*( TT(0.114713313798970)* u[j] - TT(0.290087484386815)*u[j+i] - TT(0.306681191361148)*u[j+2i] + TT(0.520262285050482)*u[j+4i]  - TT(0.051642265516119)*u[j+5i] + TT(0.013435342414630)*u[j+6i] )/Δx
    uₓ[j+4i]    = α*uₓ[j+4i]+ TT(i)*( TT(-0.036210680656541)*u[j] + TT(0.105400944933782)*u[j+i] + TT(0.015764336127392)*u[j+2i] + TT(-0.707905442575989)*u[j+3i] + TT(0.769199413962647)*u[j+5i] - TT(0.164529643265203)*u[j+6i] + T(0.018281071473911)*u[j+7i] )/Δx
    uₓ[j+5i]    = α*uₓ[j+5i]+ TT(i)*( TT(-0.011398193015050)*u[j] + TT(0.020437334208704)*u[j+i] + TT(0.011220896474665)*u[j+2i] + TT( 0.063183694641876)*u[j+3i] - TT(0.691649024426814)*u[j+4i] + TT(0.739709139060752)*u[j+6i] + TT(-0.147941827812150)*u[j+7i] + TT(0.016437980868017)*u[j+8i] )/Δx
end




@inline function FirstDerivativeBoundary!(uₓ::AT,K::AT,u::AT,Δx::TT,NT::NodeType,::Val{2},α::TT) where {TT,AT<:AbstractVector{TT}}
    NT == Left ? i = 1 : i = -1
    NT == Left ? j = 1 : j = length(u)
    uₓ[j] = α*uₓ[j] + K[j]*TT(i)*(u[j+i] - u[j])/Δx
end
@inline function FirstDerivativeBoundary!(uₓ::AT,K::AT,u::AT,Δx::TT,NT::NodeType,::Val{4},α::TT) where {TT,AT<:AbstractVector{TT}}
    NT == Left ? i = 1 : i = -1
    NT == Left ? j = 1 : j = length(u)
    uₓ[j]       = α*uₓ[j]   + K[j]      *TT(i)*(TT(-24/17)*u[j] + TT(59/34)*u[j+i]   + TT(-4/17)*u[j+2i]   + TT(-3/34)*u[j+3i])/Δx
    uₓ[j+i]     = α*uₓ[j+i] + K[j+i]    *TT(i)*(TT(-1/2)*u[j]   + TT(1/2)*u[j+2i])/Δx
    uₓ[j+2i]    = α*uₓ[j+2i]+ K[j+2i]   *TT(i)*(TT(4/43)*u[j]  + TT(-59/86)*u[j+i]   + TT(59/86)*u[j+3i]  + TT(-4/43)*u[j+4i])/Δx
    uₓ[j+3i]    = α*uₓ[j+3i]+ K[j+3i]   *TT(i)*(TT(3/98)*u[j]  + TT(-59/98)*u[j+2i]  + TT(32/49)*u[j+4i]  + TT(-4/49)*u[j+5i])/Δx
end
@inline function FirstDerivativeBoundary!(uₓ::AT,K::AT,u::AT,Δx::TT,NT::NodeType,::Val{6},α::TT) where {TT,AT<:AbstractVector{TT}}
    NT == Left ? i = 1 : i = -1
    NT == Left ? j = 1 : j = length(u)
    uₓ[j]       = α*uₓ[j]   + K[j]      *TT(i)*( TT(-1.582533518939116)*u[j] + TT(2.033378678700676)*u[j+i] - TT(0.141512858744873)*u[j+2i] + TT(-0.450398306578272)*u[j+3i] + TT(0.104488069284042)*u[j+4i] + TT(0.036577936277544)*u[j+5i] )/Δx
    uₓ[j+i]     = α*uₓ[j+i] + K[j+i]    *TT(i)*( TT(-0.462059195631158)*u[j] + TT(0.287258622978251)*u[j+2i] +TT(0.258816087376832)*u[j+3i] + TT(-0.069112065532624)*u[j+4i] - TT(0.014903449191300)*u[j+5i] )/Δx
    uₓ[j+2i]    = α*uₓ[j+2i]+ K[j+2i]   *TT(i)*( TT(0.071247104721830)* u[j] - TT(0.636451095137907)*u[j+i] + TT(0.606235523609147)*u[j+3i] + TT(-0.022902190275815)*u[j+4i] - TT(0.018129342917256)*u[j+5i] )/Δx
    uₓ[j+3i]    = α*uₓ[j+3i]+ K[j+3i]   *TT(i)*( TT(0.114713313798970)* u[j] - TT(0.290087484386815)*u[j+i] - TT(0.306681191361148)*u[j+2i] + TT(0.520262285050482)*u[j+4i]  - TT(0.051642265516119)*u[j+5i] + TT(0.013435342414630)*u[j+6i] )/Δx
    uₓ[j+4i]    = α*uₓ[j+4i]+ K[j+4i]   *TT(i)*( TT(-0.036210680656541)*u[j] + TT(0.105400944933782)*u[j+i] + TT(0.015764336127392)*u[j+2i] + TT(-0.707905442575989)*u[j+3i] + TT(0.769199413962647)*u[j+5i] - TT(0.164529643265203)*u[j+6i] + T(0.018281071473911)*u[j+7i] )/Δx
    uₓ[j+5i]    = α*uₓ[j+5i]+ K[j+5i]   *TT(i)*( TT(-0.011398193015050)*u[j] + TT(0.020437334208704)*u[j+i] + TT(0.011220896474665)*u[j+2i] + TT( 0.063183694641876)*u[j+3i] - TT(0.691649024426814)*u[j+4i] + TT(0.739709139060752)*u[j+6i] + TT(-0.147941827812150)*u[j+7i] + TT(0.016437980868017)*u[j+8i] )/Δx
end


"""
    FirstDerivativeTransposeInternal
Single node 1D first derivative transpose function.
"""
function FirstDerivativeTransposeBoundary! end
function FirstDerivativeTransposeBoundary!(dest::VT,u::AT,Δx::TT,NT::NodeType{TN},::Val{2},α::TT) where {TT,VT,AT,TN}
    TN == :Left ? i = 1 : i = -1
    TN == :Left ? j = 1 : j = length(u)

    dest[j]     = α*dest[j] + TT(i)*(TT(-1)*u[j] + TT(-1)*u[j+i])

    dest
end
function  FirstDerivativeTransposeBoundary!(dest::VT,u::AT,Δx::TT,NT::NodeType{TN},::Val{4},α::TT) where {TT,VT,AT,TN}
    TN == :Left ? i = 1 : i = -1
    TN == :Left ? j = 1 : j = length(u)
    
    dest[j]     = α*dest[j] +    TT(i)*( TT(-24//17)*u[j]+  TT(-1//2)*u[j+i] +  TT(4//43)*u[j+2i] + TT(3//98)*u[j+3i] )/Δx
    dest[j+i]   = α*dest[j+i] +  TT(i)*( TT(59//34)*u[j] +                                          TT(-59//86)*u[j+2i] )/Δx
    dest[j+2i]  = α*dest[j+2i] + TT(i)*( TT(-4//17)*u[j] +  TT(1//2)* u[j+i] +                      TT(-59//98)*u[j+3i] +   TT(-1//12)*u[j+4i] )/Δx
    dest[j+3i]  = α*dest[j+3i] + TT(i)*( TT(-3//37)*u[j] +                      TT(59//86)*u[j+2i]+                         TT(2//3)*u[j+4i] +  TT(-1//12)*u[j+5i])/Δx 
    dest[j+4i]  = α*dest[j+4i] + TT(i)*(                                        TT(-4//43)*u[j+2i]+ TT(32//49)*u[j+3i] +                       TT(2//3)*u[j+5i] +   TT(-1//12)*u[j+6i] )/Δx
    dest[j+5i]  = α*dest[j+5i] + TT(i)*(                                                            TT(-4//49)*u[j+3i] +    TT(2//3)*u[j+4i] +                      TT(2//3)*u[j+6i] + TT(-1//12)*u[j+7i] )/Δx 

    dest
end
function FirstDerivativeTransposeBoundary!(dest::VT,u::AT,c::AT,Δx::TT,NT::NodeType{TN},::Val{2},α::TT) where {TT,VT,AT,TN}
    TN == :Left ? i = 1 : i = -1
    TN == :Left ? j = 1 : j = length(u)

    dest[j]     = α*dest[j] + TT(i)*c[j]*(TT(-1)*u[j] + TT(-1)*u[j+i])

    dest
end
function  FirstDerivativeTransposeBoundary!(dest::VT,u::AT,c::AT,Δx::TT,NT::NodeType{TN},::Val{4},α::TT) where {TT,VT,AT,TN}
    TN == :Left ? i = 1 : i = -1
    TN == :Left ? j = 1 : j = length(u)
    
    dest[j]     = α*dest[j] +    TT(i)*c[j]     *( TT(-24//17)*u[j]+  TT(-1//2)*u[j+i] +  TT(4//43)*u[j+2i] + TT(3//98)*u[j+3i] )/Δx
    dest[j+i]   = α*dest[j+i] +  TT(i)*c[j+i]   *( TT(59//34)*u[j] +                                          TT(-59//86)*u[j+2i] )/Δx
    dest[j+2i]  = α*dest[j+2i] + TT(i)*c[j+2i]  *( TT(-4//17)*u[j] +  TT(1//2)* u[j+i] +                      TT(-59//98)*u[j+3i] +   TT(-1//12)*u[j+4i] )/Δx
    dest[j+3i]  = α*dest[j+3i] + TT(i)*c[j+3i]  *( TT(-3//37)*u[j] +                      TT(59//86)*u[j+2i]+                         TT(2//3)*u[j+4i] +  TT(-1//12)*u[j+5i])/Δx 
    dest[j+4i]  = α*dest[j+4i] + TT(i)*c[j+4i]  *(                                        TT(-4//43)*u[j+2i]+ TT(32//49)*u[j+3i] +                       TT(2//3)*u[j+5i] +   TT(-1//12)*u[j+6i] )/Δx
    dest[j+5i]  = α*dest[j+5i] + TT(i)*c[j+5i]  *(                                                            TT(-4//49)*u[j+3i] +    TT(2//3)*u[j+4i] +                      TT(2//3)*u[j+6i] + TT(-1//12)*u[j+7i] )/Δx 

    dest
end


function FirstDerivativeTranspose!(dest::VT,u::AT,n::Int,Δx::TT,order::Int,α::TT) where {TT,VT<:AbstractVector{TT},AT<:AbstractVector{TT}}
    order == 2 ? m = 2 : m = 7
    for i = m:n-m+1
        @inbounds dest[i] = α*dest[i] + FirstDerivativeInternal(u,Δx,Val(order),i,TT(1))
    end
    FirstDerivativeTransposeBoundary!(dest,u,Δx,Left,   Val(order),α)
    FirstDerivativeTransposeBoundary!(dest,u,Δx,Right,  Val(order),α)
    dest
end

function FirstDerivativeTranspose!(dest::VT,u::AT,c::AT,n::Int,Δx::TT,order::Int,α::TT) where {TT,VT<:AbstractVector{TT},AT<:AbstractVector{TT}}
    order == 2 ? m = 2 : m = 7
    for i = m:n-m+1
        @inbounds dest[i] = α*dest[i] + c[i]*FirstDerivativeInternal(u,Δx,Val(order),i,TT(1))
    end
    FirstDerivativeTransposeBoundary!(dest,u,c,Δx,Left,   Val(order),α)
    FirstDerivativeTransposeBoundary!(dest,u,c,Δx,Right,  Val(order),α)
    dest
end