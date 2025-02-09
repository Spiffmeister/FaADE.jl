


struct MagneticField{BT<:Union{Function,Nothing},STATE,TT<:Real,AT<:Matrix{TT}}
    B               :: BT
    istimedependant :: Bool
    Bp              :: AT
end

struct FieldLineIntercept{F<:Union{Function,Nothing}}
    Intercept :: F
end

"""
    ParallelData{TT<:Real,
        DIM,
        PGT,
        GT,
        BT,
        IT} <: ParallelGridType
"""
struct ParallelData{TT<:Real,
        DIM,
        PGT,
        GT,
        BT,
        IT,
        FLIT} <: ParallelGridType
    PGrid       :: PGT # ParallelGrid{TT,DIM,Matrix{TT}}
    κ           :: TT
    τ           :: TT
    Intercept   :: FLIT
    Interpolant :: IT
    gridx       :: GT
    gridy       :: GT
    Δx          :: TT
    Δy          :: TT

    H           :: CompositeH{2,TT,Vector{TT},DiagonalH{TT,Vector{TT}}}

    w_f         :: Matrix{TT}
    w_b         :: Matrix{TT}

    u           :: Matrix{TT}

    MagneticField   :: BT

    τ_i         :: Vector{TT}
    
end
"""
    ParallelData(PGrid::ParallelGrid,G::Grid1D{TT},order::Int;κ=TT(1),intercept=nothing) where TT
"""
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

    u = zeros(TT,1,1)

    Bp = zeros(TT,(1,1))
    MF = MagneticField{Nothing,:EQUILIBRIUM,TT,typeof(Bp)}(nothing,false,Bp)
    
    ParallelData{TT,1,typeof(PGrid),typeof(gridx),typeof(MF),Nothing}(PGrid,κ,τ,intercept_fieldlines,nothing,gridx,gridy,Δx,Δy,H,w_f,w_b,u,MF,[TT(0)])
end
"""
    ParallelData(PGrid::ParallelGrid,G::Grid2D{TT},order::Int;κ=TT(1),intercept=nothing,B=nothing,interpolant=nothing,remap=nothing) where {TT}
"""
function ParallelData(PGrid::ParallelGrid,G::Grid2D{TT},order::Int;κ::TT=TT(1),intercept=nothing,interpolant=nothing,periodicy=false,B=nothing) where {TT}

    τ = TT(1)

    Hx = DiagonalH(order,G.Δx,G.nx)
    Hy = DiagonalH(order,G.Δy,G.ny)

    H = CompositeH(Hx,Hy)

    # if isnothing(remap)
    gridx = G.gridx[1:end,1]
    gridy = G.gridy[1,1:end]
    # else
    #     gridx = remap(G.gridx,G.gridy)[1:end,1]
    #     gridy = remap(G.gridx,G.gridy)[1,1:end]
    # end

    w_f = zeros(TT,size(G))
    w_b = zeros(TT,size(G))

    u = zeros(TT,1,1)

    Bp = zeros(TT,(3,3))
    MF = MagneticField{typeof(B),:EQUILIBRIUM,TT,typeof(Bp)}(B,false,Bp)

    if isnothing(B)
        # @warn "B not provided, perpendicular solve may not be performed correctly."
    end

    if typeof(interpolant) <: Symbol
    #     if interpolant == :bicubic
    #         Interpolator(u) = BicubicInterpolator(G.gridx,G.gridy,u)
    #     elseif interpolant == :bilinear
    #         Interpolator(u) = LinearInterpolator(G.gridx,G.gridy,u)
    #     end
        if interpolant == :chs

            nx = G.nx
            ny = G.ny

            if periodicy
                # gridx = TT(0):G.Δx:TT(1)
                # gridy = TT(0):G.Δy:TT(1)
                gridx = G.gridx[:,1:ny-1][:]
                gridy = G.gridy[:,1:ny-1][:]
            else
                gridx = G.gridx[:]
                gridy = G.gridy[:]
            end

            z = zeros(TT,length(gridx))
            dzdx = zeros(TT,length(gridx))
            dzdy = zeros(TT,length(gridx))
            Interpolator = BivariateCHSInterpolation(gridx,gridy,z,dzdx,dzdy)

            if isa(intercept,Function)
                Intercept = intercept
            elseif isnothing(intercept)
                Intercept = nothing
            end


            gridx = G.gridx[1:end,1]
            gridy = G.gridy[1,1:end]
        end
    elseif isnothing(interpolant)
        Interpolator = interpolant
    end

    ParallelData{TT,2,typeof(PGrid),typeof(gridx),typeof(MF),typeof(Interpolator),typeof(Intercept)}(PGrid,κ,τ,Intercept,Interpolator,gridx,gridy,G.Δx,G.Δy,H,w_f,w_b,u,MF,[TT(0.0)])
end



"""
    ParallelMultiBlock{TT,DIM,IT}
Stores the parallel data for multiblock problems.
"""
struct ParallelMultiBlock{TT<:Real,
        DIM,
        IT,
        CT} <: ParallelGridType
    
    PData       :: Dict{Int,ParallelData{TT,DIM}}
    uglobal     :: Vector{Matrix{TT}}
    Interpolant :: IT
    Intercept   :: CT
    τ           :: Vector{TT}
end
function ParallelMultiBlock(PGrid::Dict,G::GridMultiBlock{TT},order::Int;κ=TT(1),interpopts=Dict()) where {TT}

    intercept = nothing
    interpolation_mode = nothing
    periodicy = false
    B = nothing

    if haskey(interpopts,"intercept")
        intercept = interpopts["intercept"]
    end
    if haskey(interpopts,"interpolant")
        interpolation_mode = interpopts["interpolant"]
    end
    if haskey(interpopts,"periodicy")
        periodicy = interpopts["periodicy"]
    end
    if haskey(interpopts,"B")
        B = interpopts["B"]
    end


    PData  = Dict{Int,ParallelData}()
    for I in eachgrid(G)
        PData[I] = ParallelData(PGrid[I],G.Grids[I],order,κ=κ,intercept=intercept,interpolant=interpolation_mode,periodicy=periodicy,B=B)
    end
    # ParallelData = Tuple(collect(values(PData)))

    # interpolant = [(u->LinearInterpolator(G.Grids[I].gridx,G.Grids[I].gridy,u)) for I in eachgrid(G)]
    # interpolant = :nearest
    # if interpolant == :CHS
    # @show PData[3].Interpolant.x[1:10]
    # @show PData[3].Interpolant.y[1:10]
    # TmpInterpolant = [PData[I].Interpolant for I in eachindex(PData)]
    # TmpInterpolant = [PData[I].Interpolant for I in keys(PData)]
    # @show PData[1].Interpolant
    
    TmpInterpolant = []
    for I in 1:length(PData)
        push!(TmpInterpolant,PData[I].Interpolant)
        # @show PData[I].Interpolant
    end
    # @show TmpInterpolant[1]

    # TmpInterpolant = [PData[I].Interpolant for I in keys()]
    Interp = Tuple(TmpInterpolant)
    # end
    if isa(intercept,Function)
        TmpIntercept = [PData[I].Intercept for I in eachindex(PData)]
        intercept = Tuple(TmpIntercept)
    end

    uglobal = [zeros(TT,size(G.Grids[I])) for I in eachgrid(G)]

    # u = zeros(TT,size(GridMultiBlock))

    τglobal = zeros(TT,length(Interp))

    return ParallelMultiBlock{TT,2,typeof(Interp),typeof(intercept)}(PData,uglobal,Interp,intercept,τglobal)
end




"""
"""
function construct_parallel(PGrid::ParallelGrid,grid::Grid2D{TT},order::Int,interpoptions::Dict=Dict()) where TT


    keyopts = keys(interpoptions)

    if "coefficient" ∈ keyopts
        κ = interpoptions["coefficient"]
    else
        κ = TT(1)
    end

    if "field" ∈ keyopts
        B = interpoptions["field"]
    else
        B = nothing
    end

    if "periodicy" ∈ keyopts
        periodicy = interpoptions["periodicy"]
    else
        periodicy = false
    end

    if "interpolant" ∈ keyopts
        interpolant = interpoptions["interpolant"]
        if (interpolant == :chs) && !("boundaryconditions" ∈ keyopts) && !(("bounds" ∈ keyopts) || ("intercept" ∈ keyopts))
            error("If CubicHermiteSpline interpolation is used the boundary conditions and outer domain boundaries must also be supplied by the \"boundaryconditions\" and \"bounds\" OR \"intercept\" Dict keys.")
        else
            intercept = interpoptions["boundaryconditions"]

            if !periodicy #double  check if y is periodic or not since this can break if it should be
                @warn "CubicHermiteSpline interpolation specified, if the domain is periodic then you should specify `\"periodicy\"=>true`."
            end
        end
    else
        interpolant = :none
        intercept = nothing
    end

    

    PData = ParallelData(PGrid,grid,order,κ=κ,B=B,intercept=intercept,interpolant=interpolant,periodicy=periodicy)
    return PData
end