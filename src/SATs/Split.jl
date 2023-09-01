


"""
    SAT_Split

TODO: Write functinos and struct for split domain SATs
"""
struct SAT_Split <: SimultanousApproximationTerm{:Interface} end




"""
Split_domain(u⁻::Vector{Float64},u⁺::Vector{Float64},Δx⁻::Float64,Δx⁺::Float64,c⁻,c⁺;
    order::Int64=2,order⁻::Int64=2,order⁺::Int64=2,separate_forcing::Bool=false)
Simulatenous approximation term for a split domain
TODO: Fix for when `order⁺ != order⁻`
"""
function Split_domain(u⁻::Vector{Float64},u⁺::Vector{Float64},Δx⁻::Float64,Δx⁺::Float64,c⁻,c⁺;
        order::Int64=2,order⁻::Int64=2,order⁺::Int64=2,separate_forcing::Bool=false)

    h⁻ = hval(order⁻)
    h⁺ = hval(order⁺)

    # Left side
    α₀ = -0.5
    τ₁ = 0.5
    τ₀ = max(c⁻[end]/(2*h⁻*Δx⁻), c⁺[1]/(2*h⁺*Δx⁺))

    SAT₀ = zeros(Float64,order⁻+order⁺)
    SAT₁ = zeros(Float64,order⁻+order⁺)
    SAT₂ = zeros(Float64,order⁻+order⁺)

    F = zeros(Float64,order⁻+order⁺)

    # Function condition
    L₀u = zeros(Float64,order⁻+order⁺)
    L₀u[order⁻] += u⁻[end]
    L₀u[order⁻+1] += u⁻[end]

    SAT₀[1:order⁻] = τ₀/(h⁻ * Δx⁻) * L₀u[1:order⁻]
    SAT₀[order⁻+1:end] = τ₀/(h⁺ * Δx⁺) * L₀u[order⁻+1:end]

    F[1:order⁻] = -τ₀/(h⁻ * Δx⁻) * u⁺[1]
    F[order⁻+1:end] = τ₀/(h⁺ * Δx⁺) * u⁻[end]


    # Symmeteriser
    D₁ᵀL₀u⁻ = boundary_D₁ᵀ(L₀u[1:order⁻],Δx⁻,order⁻)
    D₁ᵀL₀u⁺ = boundary_D₁ᵀ(L₀u[order⁻+1:end],Δx⁺,order⁺)

    SAT₁[1:order⁻] = τ₁/(h⁻ * Δx⁻) * c⁻[end] * D₁ᵀL₀u⁻[order⁻+1:end]
    SAT₁[order⁻+1:end] = -τ₁/(h⁺ * Δx⁺) * c⁺[1] * D₁ᵀL₀u⁺[1:order⁺]

    # Derivative condition
    D₁u⁻ = boundary_D₁(u⁻,Δx⁻,order⁻)
    D₁u⁺ = boundary_D₁(u⁺,Δx⁺,order⁺)

    SAT₂[1:order⁻] = α₀/(h⁻ * Δx⁻) * c⁻[end] * D₁u⁻[order⁻+1:end]
    SAT₂[order⁻+1:end] = α₀/(h⁺ * Δx⁺) * c⁺[1] * D₁u⁺[1:order⁺]

    SAT = SAT₀ + SAT₁ + SAT₂

    SAT⁻ = SAT[1:order⁻]
    SAT⁺ = SAT[order⁻+1:end]


    if !separate_forcing
        SAT⁻ += F[1:order⁻]
        SAT⁺ += F[order⁻+1:end]
        return SAT⁻, SAT⁺
    else
        return SAT⁻, SAT⁺, F[1:order⁻], F[order⁻+1:end]
    end

end



struct SAT_Interface{
        TN<:NodeType,
        TT<:Real,
        TV<:Vector{TT},
        F1<:Function} <: SimultanousApproximationTerm{:Interface}

    side    :: TN
    axis    :: Int
    order   :: Int
    D₁ᵀE₀   :: TV
    D₁ᵀEₙ   :: TV
    D₁E₀    :: TV
    D₁Eₙ    :: TV
    τ₀  :: TT
    α₀  :: TT
    τ₀  :: TT
    loopaxis :: F1

    function SAT_Split(Δx::TT,side::TN,axis::Int,order::Int) where {TT,TN}
        D₁ᵀE₀ = BoundaryDerivativeTranspose(Left,order,Δx)
        D₁ᵀEₙ = BoundaryDerivativeTranspose(Right,order,Δx)
        E₀D₁ = BoundaryDerivative(Left,Δx,order)
        EₙD₁ = BoundaryDerivative(Right,Δx,order)

        α₀, τ₁, τ₀ = SATpenalties(Interface,Δx,order)

        loopaxis = SelectLoopDirection(axis)

        new{TN,TT,Vector{TT},typeof(loopaxis)}(side,axis,order,
            D₁ᵀE₀,D₁ᵀEₙ,E₀D₁,EₙD₁,τ₀,α₀,τ₁,loopaxis)
    end

end


function SAT_Interface!(u⁻::Vector{Float64},u⁺::Vector{Float64},Δx⁻::Float64,Δx⁺::Float64,c⁻,c⁺;
        order::Int64=2,order⁻::Int64=2,order⁺::Int64=2,separate_forcing::Bool=false)

    for (S⁻,S⁺,U⁻,U⁺,K⁻,K⁺) in zip(loopaxis(u),)

        
        S⁻[end] = S⁻[end]   + τ₀/(h⁻ * Δx⁻) * U⁻[end] #SAT₀
        S⁺[1]   = S⁺[1]     + τ₀/(h⁺ * Δx⁺) * U⁺[1] #SAT₀
        
        for i = 1:order⁻
            S⁻[i] = S⁻[i] + τ₁/(h⁻ * Δx⁻) * K⁻[end] * D₁ᵀE₀[i] * U⁻[end-order⁻+i]
            S⁻[i] = S⁻[i] + α₀/(h⁻ * Δx⁻) * K⁻[end] * D₁E₀[i] * U⁻[end-order⁻+i]

            S⁺[i] = S⁺[i] + τ₁/(h⁺ * Δx⁺) * K⁺[1] * D₁ᵀE₀[i] * U⁺[i]
            S⁺[i] = S⁺[i] + α₀/(h⁺ * Δx⁺) * K⁺[1] * D₁E₀[i] * U⁺[i]
        end

    end

end



function (SI::SAT_Interface{NodeType{:Left,DIM},TT})(cache::AT,c::AT,rhs::AT) where {TT,DIM,AT}
    # rhs = u⁻ - u⁺
    for (S⁺,U⁻,K⁺) in zip(SI.loopaxis(cache),SI.loopaxis(rhs),SI.loopaxis(c))
        for i = 1:SI.order
            # τ₁ K D₁ᵀL₀ u + αL₀KD₁u
            S⁺[i] += K⁺[1] * (SI.τ₁*SI.D₁ᵀE₀[i] + SI.α₀*SI.D₁E₀[i]) * (U⁻[i] - S⁺[i])
        end
        S⁺[1] += -τ₀ * (U⁻[1] - S⁺[1]) # L₀u = u⁻ - u⁺
    end
end
function (SI::SAT_Interface{NodeType{:Right,DIM},TT})(cache::AT,c::AT,rhs::AT) where {TT,DIM,AT}
    for (S⁻,U⁺,K⁻) in zip(SI.loopaxis(cache),SI.loopaxis(rhs),SI.loopaxis(c))
        for i = 1:SI.order
            S⁻[i] += SI.τ₁ * K⁻[end] * (SI.D₁ᵀE₀[i] + SI.D₁E₀[i]) * U⁻[end-SI.order+i]
        end
        S⁻[end] += SI.τ₀ * U⁻[end]
    end
end
