



"""
    GetAxis
Return the axis along which the node lies
"""
GetAxis(::NodeType{T,Ax}) where {T,Ax} = Int(Ax)

"""
    GetDim
Return the data structures associated dimension
"""
function GetDim end
GetDim(::BoundaryStorage{T,D}) where {T,D} = Int(D)
GetDim(::DataBlock{T,D}) where {T,D} = Int(D)
GetDim(::GridType{T,D}) where {T,D} = Int(D)

"""
    GetMinΔ
Return the miniumum
"""
function GetMinΔ end
GetMinΔ(grid::Grid1D) = grid.Δx
GetMinΔ(grid::Grid2D) = min(grid.Δx,grid.Δy)
