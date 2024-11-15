
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

