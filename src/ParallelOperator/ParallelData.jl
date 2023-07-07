
"""
    ParGrid{TT,AT<:AbstractArray{TT}}
Storing the x-y coordinates of a parallel grid
"""
struct ParGrid{TT,AT<:AbstractArray{TT}}
    x   :: AT
    y   :: AT
end

"""
    ParallelGrid{TT,DIM,AT,PA}
Stores the current, forward and backward planes for the parallel tracing.

In 2D arrays are of format ``[(x_1,y_1),(x_1,y_2),...,(x_n,y_n)]``
"""
struct ParallelGrid{TT<:Real,
        DIM,
        AT<:Matrix{TT},
        PA<:Matrix{Vector{TT}}}

    plane   :: PA
    Bplane  :: ParGrid{TT,AT}
    Fplane  :: ParGrid{TT,AT}
    # w_f     :: AT
    # w_b     :: AT
end



struct FieldLineIntercept{F<:Union{Function,Nothing}}
    Intercept :: F
end


struct ParallelData{TT<:Real,
        DIM,
        FLG}
    PGrid       :: ParallelGrid{TT,DIM,Matrix{TT},Matrix{Vector{TT}}}
    κ           :: TT
    τ           :: TT
    Intercept   :: FieldLineIntercept
    gridx       :: LinRange{TT,Int64}
    gridy       :: LinRange{TT,Int64}
    Δx          :: TT
    Δy          :: TT

    function ParallelData(PGrid::ParallelGrid,Grid::Grid2D{TT};κ=TT(1),intercept=nothing,fieldlinegeom=:confined) where {TT}

        intercept_fieldlines = FieldLineIntercept(intercept)
        # intercept_fieldlines = intercept

        # K = DiffusionCoefficient(κ)

        τ = -sqrt((P.gridx[end]-P.gridx[1])*(P.gridy[end]-P.gridy[1])/(P.Δx*P.Δy))


        gridx = LinRange(Grid.gridx[1],Grid.gridx[end],Grid.nx)
        gridy = LinRange(Grid.gridy[1],Grid.gridy[end],Grid.ny)

        new{TT,2,Mode,fieldlinegeom}(PGrid,κ,τ,intercept_fieldlines,gridx,gridy,Grid.Δx,Grid.Δy)
    end
end


