




"""
    Grid1D
Grid data structure for 1 dimensional problems.
    
    `Grid1D{T}(𝒟::Vector,n::Integer)`
Inputs:
- Vector of domain boundaries `[x_0,x_N]`
- Number of grid points.

Returns:
- Struct for 1D grid object containing vector of grid points, ``\\Delta x`` and ``n``.
"""
struct Grid1D{T} <: GridType{T,1}
    grid    :: Vector{T}
    Δx      :: T
    n       :: Integer
    function Grid1D(𝒟::Vector{T},n::Integer) where T
        Δx = (𝒟[2]-𝒟[1])/(n-1)
        x = collect(range(𝒟[1],𝒟[2],length=n))
        new{T}(x,Δx,n)
    end
end



"""
    Grid2D
Grid data structure for 2 dimensional problems.

    `Grid2D{T}(𝒟x::Vector,𝒟y::Vector,nx::Integer,ny::Integer)`
Inputs:
- Domain boundaries in ``x``
- Domain boundaries in ``y``
- Number of nodes in ``x``
- Number of nodes in ``y``

Returns:
- Struct for 2D grid object containing grid points in ``x`` and ``y``, ``\\Delta x`` and ``\\Delta y``, and ``n_x`` and ``n_y``.
"""
struct Grid2D{T} <: GridType{T,2}
    gridx   :: Vector{T}
    gridy   :: Vector{T}
    Δx      :: T
    Δy      :: T
    nx      :: Integer
    ny      :: Integer
    function Grid2D(𝒟x::Vector{T},𝒟y::Vector{T},nx::Integer,ny::Integer) where T
        gx = Grid1D(𝒟x,nx)
        gy = Grid1D(𝒟y,ny)
    
        new{T}(gx.grid,gy.grid, gx.Δx,gy.Δx, gx.n,gy.n)
    end
end






"""
    GetMinΔ
Return the miniumum grid size between ``\\Delta x`` and ``\\Delta y``
"""
function GetMinΔ end
GetMinΔ(grid::Grid1D) = grid.Δx
GetMinΔ(grid::Grid2D) = min(grid.Δx,grid.Δy)
