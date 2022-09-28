


abstract type GridType{T<:AbstractFloat,N} end


"""
    Grid1D{T}
Grid data structure for 1 dimensional problems
"""
struct Grid1D{T} <: GridType{T,1}
    grid    :: Vector{T}
    Δx      :: T
    n       :: Int
    function Grid1D(𝒟::Vector{T},n::Int) where T
        Δx = (𝒟[2]-𝒟[1])/(n-1)
        x = collect(range(𝒟[1],𝒟[2],length=n))
        new{T}(x,Δx,n)
    end
end



"""
    Grid2D{T}
Grid data structure for 2 dimensional problems
"""
struct Grid2D{T} <: GridType{T,2}
    gridx   :: AbstractArray{T}
    gridy   :: AbstractArray{T}
    Δx      :: T
    Δy      :: T
    nx      :: Int
    ny      :: Int
    function Grid2D(𝒟x::Vector{T},𝒟y::Vector{T},nx::Int,ny::Int) where T
        gx = Grid1D(𝒟x,nx)
        gy = Grid1D(𝒟y,ny)
    
        new{T}(gx.grid,gy.grid, gx.Δx,gy.Δx, gx.n,gy.n)
    end
end





