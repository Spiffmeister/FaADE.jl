


abstract type DerivativeOperatorType{DIM} end

struct DerivativeOperator1D{Diffuse,Advect} <: DerivativeOperatorType{1} end
struct DerivativeOperator2D{DiffuseX,DiffuseY,AdvectX,AdvectY} <: DerivativeOperatorType{2} end


struct DerivativeOrder{O} end
struct DerivativeOperator{OrderX,OrderY,MixedOrderX,MixedOrderY} <: DerivativeOperatorType{2} end



# Base.show(io::IO,DO::DerivativeOrder{O}) where O = print(io, "order ",O," second derivative operator")
