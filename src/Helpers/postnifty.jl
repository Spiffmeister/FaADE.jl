



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
