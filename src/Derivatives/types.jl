


abstract type DerivativeOperatorType{DIM} end

struct DerivativeOperator1D{Diffuse,Advect} <: DerivativeOperatorType{1} end
struct DerivativeOperator2D{DiffuseX,DiffuseY,AdvectX,AdvectY} <: DerivativeOperatorType{2} end


struct DerivativeOrder{O} end
struct DerivativeOperator{TT<:Real,
        DIM,
        DO<:DerivativeOrder,
        COEFF} <: DerivativeOperatorType{DIM}
    order   :: DO
    nx      :: Int64
    ny      :: Int64
    Δx      :: TT
    Δy      :: TT
end

GetOrder(D::DerivativeOrder{O}) where {O} = O
GetOrder(O::Int) = O

Base.show(io::IO,DO::DerivativeOperator{TT,DIM,O}) where {TT,DIM,O} = print("Order ",GetOrder(DO.order)," ",DIM," dimensional diffusion SBP operator.")


# Base.show(io::IO,DO::DerivativeOrder{O}) where O = print(io, "order ",O," second derivative operator")
