
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
        AT<:Matrix{TT}}

    Bplane  :: ParGrid{TT,AT}
    Fplane  :: ParGrid{TT,AT}
end



struct FieldLineIntercept{F<:Union{Function,Nothing}}
    Intercept :: F
end


struct ParallelData{TT<:Real,
        DIM,
        GT}
    PGrid       :: ParallelGrid{TT,DIM,Matrix{TT}}
    κ           :: TT
    τ           :: TT
    Intercept   :: FieldLineIntercept
    gridx       :: GT
    gridy       :: GT
    Δx          :: TT
    Δy          :: TT

    H           :: CompositeH{2,TT,Vector{TT},:Diagonal}

    w_f         :: Matrix{TT}
    w_b         :: Matrix{TT}

    function ParallelData(PGrid::ParallelGrid,G::Grid2D{TT};κ=TT(1),intercept=nothing,order=2) where {TT}

        intercept_fieldlines = FieldLineIntercept(intercept)
        # intercept_fieldlines = intercept

        # K = DiffusionCoefficient(κ)

        # τ = -sqrt(1.0/((G.gridx[end]-G.gridx[1])*(G.gridy[end]-G.gridy[1]) * G.Δx*G.Δy))
        # τ = -1/(G.Δx*G.Δy)
        # τ = TT(-1)

        τ = TT(1)

        Hx = DiagonalH(order,G.Δx,G.nx)
        Hy = DiagonalH(order,G.Δy,G.ny)

        H = CompositeH(Hx,Hy)


        gridx = LinRange(G.gridx[1],G.gridx[end],G.nx)
        gridy = LinRange(G.gridy[1],G.gridy[end],G.ny)

        w_f = zeros(TT,size(G))
        w_b = zeros(TT,size(G))

        new{TT,2,typeof(gridx)}(PGrid,κ,τ,intercept_fieldlines,gridx,gridy,G.Δx,G.Δy,H,w_f,w_b)
    end
end


