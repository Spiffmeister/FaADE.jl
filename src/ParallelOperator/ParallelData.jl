
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
        GT}
    PGrid       :: ParallelGrid{TT,DIM,Matrix{TT},Matrix{Vector{TT}}}
    κ           :: TT
    τ           :: TT
    Intercept   :: FieldLineIntercept
    gridx       :: GT
    gridy       :: GT
    Δx          :: TT
    Δy          :: TT

    w_f         :: Matrix{TT}
    w_b         :: Matrix{TT}

    function ParallelData(PGrid::ParallelGrid,G::Grid2D{TT};κ=TT(1),intercept=nothing) where {TT}

        intercept_fieldlines = FieldLineIntercept(intercept)
        # intercept_fieldlines = intercept

        # K = DiffusionCoefficient(κ)

        τ = -sqrt(1.0/((G.gridx[end]-G.gridx[1])*(G.gridy[end]-G.gridy[1]) * G.Δx*G.Δy))


        gridx = LinRange(G.gridx[1],G.gridx[end],G.nx)
        gridy = LinRange(G.gridy[1],G.gridy[end],G.ny)

        w_f = zeros(TT,size(G))
        w_b = zeros(TT,size(G))

        new{TT,2,typeof(gridx)}(PGrid,κ,τ,intercept_fieldlines,gridx,gridy,G.Δx,G.Δy,w_f,w_b)
    end
end


