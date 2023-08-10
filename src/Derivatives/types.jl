


abstract type DerivativeOperatorType{DIM} end

struct DerivativeOperator1D{Diffuse,Advect} <: DerivativeOperatorType{1} end
struct DerivativeOperator2D{DiffuseX,DiffuseY,AdvectX,AdvectY} <: DerivativeOperatorType{2} end


struct DerivativeOrder{O} end
struct DerivativeOperator{TT,DIM,Diffuse,Mixed,Advect} <: DerivativeOperatorType{DIM} 
    order   :: DerivativeOrder
    nx      :: Int64
    ny      :: Int64
    Δx      :: TT
    Δy      :: TT
end

GetOrder(D::DerivativeOrder{O}) where {O} = O
GetOrder(O::Int) = O




# Base.show(io::IO,DO::DerivativeOrder{O}) where O = print(io, "order ",O," second derivative operator")
