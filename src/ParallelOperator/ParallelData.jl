

abstract type ParallelGridType{TT} end



"""
    ParGrid{TT,AT<:AbstractArray{TT}}
Storing the x-y coordinates of a parallel grid
"""
struct ParGrid{TT,AT<:AbstractArray{TT}}
    x       :: AT
    y       :: AT
    subgrid :: Matrix{Int}
end

"""
    ParallelGrid{TT,DIM,AT,PA}
Stores the current, forward and backward planes for the parallel tracing.

In 2D arrays are of format ``[(x_1,y_1),(x_1,y_2),...,(x_n,y_n)]``
"""
struct ParallelGrid{TT<:Real,
        DIM,
        AT<:Matrix{TT}}

    Bplane  :: ParGrid{TT,AT}
    Fplane  :: ParGrid{TT,AT}
end

struct MagneticField{BT<:Union{Function,Nothing},STATE,TT<:Real,AT<:Matrix{TT}}
    B               :: BT
    istimedependant :: Bool
    Bp              :: AT
end

struct FieldLineIntercept{F<:Union{Function,Nothing}}
    Intercept :: F
end


struct ParallelData{TT<:Real,
        DIM,
        GT,
        BT,
        IT}
    PGrid       :: ParallelGrid{TT,DIM,Matrix{TT}}
    κ           :: TT
    τ           :: TT
    Intercept   :: FieldLineIntercept
    Interpolant :: IT
    gridx       :: GT
    gridy       :: GT
    Δx          :: TT
    Δy          :: TT

    H           :: CompositeH{2,TT,Vector{TT},DiagonalH{TT,Vector{TT}}}

    w_f         :: Matrix{TT}
    w_b         :: Matrix{TT}

    MagneticField   :: BT

    τ_i         :: Vector{TT}
    
end
function ParallelData(PGrid::ParallelGrid,G::Grid2D{TT},order::Int;κ=TT(1),intercept=nothing,B=nothing,interpolant=nothing) where {TT}

    intercept_fieldlines = FieldLineIntercept(intercept)

    τ = TT(1)

    Hx = DiagonalH(order,G.Δx,G.nx)
    Hy = DiagonalH(order,G.Δy,G.ny)

    H = CompositeH(Hx,Hy)

    gridx = G.gridx[1:end,1]
    gridy = G.gridy[1,1:end]

    w_f = zeros(TT,size(G))
    w_b = zeros(TT,size(G))

    Bp = zeros(TT,(3,3))
    MF = MagneticField{typeof(B),:EQUILIBRIUM,TT,typeof(Bp)}(B,false,Bp)

    if isnothing(B)
        # @warn "B not provided, perpendicular solve may not be performed correctly."
    end

    if typeof(interpolant) <: Nothing
        Interpolator = nothing
    elseif typeof(interpolant) <: Symbol
        if interpolant == :bicubic
            Interpolator(u) = BicubicInterpolator(G.gridx,G.gridy,u)
        elseif interpolant == :bilinear
            Interpolator(u) = LinearInterpolator(G.gridx,G.gridy,u)
        end
    end

    ParallelData{TT,2,typeof(gridx),typeof(MF),typeof(Interpolator)}(PGrid,κ,τ,intercept_fieldlines,Interpolator,gridx,gridy,G.Δx,G.Δy,H,w_f,w_b,MF,[TT(0)])
end
function ParallelData(PGrid::ParallelGrid,G::Grid1D{TT},order::Int;κ=TT(1),intercept=nothing) where TT

    intercept_fieldlines = FieldLineIntercept(intercept)

    τ = TT(-1)

    Hx = DiagonalH(order,G.Δx,G.n)

    H = CompositeH(Hx,Hx)

    gridx = G.grid
    gridy = zeros(TT,1)

    Δx = G.Δx
    Δy = TT(0)

    w_f = zeros(TT,(length(G),1))
    w_b = zeros(TT,(length(G),1))

    Bp = zeros(TT,(1,1))
    MF = MagneticField{Nothing,:EQUILIBRIUM,TT,typeof(Bp)}(nothing,false,Bp)
    
    ParallelData{TT,1,typeof(gridx),typeof(MF),Nothing}(PGrid,κ,τ,intercept_fieldlines,nothing,gridx,gridy,Δx,Δy,H,w_f,w_b,MF,[TT(0)])
end



function ParallelMultiBlock(PGrid::Dict,G::GridMultiBlock{TT},order::Int;κ=TT(1),intercept=nothing,B=nothing,interpolant=nothing) where {TT}

    PData  = Dict{Int,ParallelData}()
    for I in eachgrid(G)
        PData[I] = ParallelData(PGrid[I],G.Grids[I],order,κ=κ,intercept=intercept,B=B,interpolant=interpolant)
    end
    return PData
end

