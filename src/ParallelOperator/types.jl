
"""
    ParallelGridType
Abstract type for storage of 
"""
abstract type ParallelGridType end

"""
    ParallelMapType
"""
abstract type ParallelMapType end

"""
    Options for constructing the parallel grid
"""
mutable struct ParallelGridOptions
    periodic_x              :: Bool
    periodic_y              :: Bool
    interpolation_mode      :: Symbol
    coordinate_map          :: Union{Nothing,Function,Tuple}
    inverse_coordinate_map  :: Union{Nothing,Function}
    boundary_x              :: Union{Nothing,Function,Vector}
    boundary_y              :: Union{Nothing,Function,Vector}
    ParallelGridOptions(;periodic_x=false,periodic_y=true,interpolation_mode=:bicubic,coordinate_map=nothing,inverse_coordinate_map=nothing,boundary_x=nothing,boundary_y=nothing) = 
        new(periodic_x,periodic_y,interpolation_mode,coordinate_map,inverse_coordinate_map,boundary_x,boundary_y)
end



mutable struct ParallelMapOptions
    Intercept       :: Union{Nothing,Function}
    Interpolant     :: Union{Symbol,Function}
    MagneticField   :: Union{Nothing,Function}
    ParallelMapOptions(;Intercept=nothing,Interpolant=:CubicHermiteSpline,MagneticField=nothing) = 
        new(Intercept,Interpolant,MagneticField)
end




"""
    ParGrid{TT,AT<:AbstractArray{TT}}
Storing the x-y coordinates of a parallel grid
"""
struct ParGrid{TT,AT<:AbstractArray{TT}} <: ParallelMapType
    x       :: AT
    y       :: AT
    subgrid :: Matrix{Int}
end
"""
    ParGridLinear{TT,AT<:AbstractArray{TT}}
Storage for a parallel map using linear interpolation.
"""
struct ParGridLinear{TT,AT<:AbstractArray{TT},METHOD} <: ParallelMapType
    weight11 :: AT
    weight12 :: AT
    weight21 :: AT
    weight22 :: AT

    i       :: Matrix{Int}
    j       :: Matrix{Int}

    subgrid :: Matrix{Int}
end

"""
    ParallelGrid{TT,DIM,AT,PA}
Stores the current, forward and backward planes for the parallel tracing.

In 2D arrays are of format ``[(x_1,y_1),(x_1,y_2),...,(x_n,y_n)]``
"""
struct ParallelGrid{TT<:Real,
        DIM,
        PMT<:ParallelMapType,
        AT<:Matrix{TT}} <: ParallelGridType

    Bplane  :: PMT
    Fplane  :: PMT
end


#
# Below are some functions for IO so the parallel grid does not have to be computed every time
#   or so they can be precomputed and stored.
#

"""
    savegrid(fname,pargrid::ParallelGridType)

Saves the `ParallelGrid` object containing the backward and forward plane information to a `jld2` file with path and filename specified by `fname`.
"""
function savegrid(fname::String,pargrid::PGT)  where {PGT<:ParallelGridType}
    if PGT <: ParallelGrid
        jldsave(fname;
            Bx = pargrid.Bplane.x,
            By = pargrid.Bplane.y,
            Bsubgrid = pargrid.Bplane.subgrid,
            Fx = pargrid.Fplane.x,
            Fy = pargrid.Fplane.y,
            Fsubgrid = pargrid.Fplane.subgrid,
            type = :ParGrid)
    end
end
"""
    readgrid(fname::String)

Read the `ParallelGrid` saved by [`savegrid`](@ref) into a `ParallelGridType` object.
"""
function readgrid(fname::String)
    jldfile = jldopen(fname)
    if jldfile["type"] == :ParGrid
        x = jldfile["Bx"]
        y = jldfile["By"]
        subgrid = jldfile["Bsubgrid"]
        Bplane = ParGrid{eltype(x),typeof(x)}(x,y,subgrid)
        x = jldfile["Fx"]
        y = jldfile["Fy"]
        subgrid = jldfile["Fsubgrid"]
        Fplane = ParGrid{eltype(x),typeof(x)}(x,y,subgrid)
        return ParallelGrid{eltype(x),ndims(x),typeof(Bplane),typeof(Bplane.x)}(Bplane,Fplane)
    end
end
function readgrids(fpath::String)
    fnames = readdir(fpath)
    gdata = Dict()
    for (I,fname) in enumerate(fnames)
        gdata[I] = readgrid(fname)
    end
    return gdata
end
