


abstract type DerivativeOperatorType{DIM} end


struct DerivativeOrder{O} end


struct DerivativeOperator{TT<:Real,
        DIM,
        DO<:DerivativeOrder,
        COEFF} <: DerivativeOperatorType{DIM}
    order       :: DO
    nx          :: Int64
    ny          :: Int64
    Δx          :: TT
    Δy          :: TT
    xperiodic   :: Bool
    yperiodic   :: Bool
end


struct DiffusionOperator{TT<:Real,DO<:DerivativeOrder,COEFF} <: DerivativeOperatorType{1}
    order   :: DO
    n       :: Int64
    Δx      :: TT
    periodic:: Bool
    function DiffusionOperator(n::Int64,Δx::TT,order::Int,periodic::Bool,coeff::Symbol) where {TT<:Real}
        DO = DerivativeOrder{order}()
        return new{TT,typeof(DO),coeff}(DO,n,Δx,periodic)
    end
end
struct DiffusionOperatorND{TT,DIM,DO,COEFF,AT<:AbstractArray{TT,DIM}} <: DerivativeOperatorType{DIM}
    DO      :: NTuple{DIM,DiffusionOperator{TT,DO,COEFF}}
    cache   :: AT
    function DiffusionOperatorND(D::DiffusionOperator{TT,DO,COEFF}...) where {TT<:Real,DO,COEFF}
        DIM = length(D)
        # cache = zeros(AT,(DO[1].n,DO[2].n))
        if COEFF == :Constant
            cache = zeros(TT,(1,1))
        elseif COEFF == :Variable
            cache = zeros(TT,(D[1].n,D[2].n))
        end
        
        new{TT,DIM,typeof(D[1].order),COEFF,typeof(cache)}(D,cache)
    end
end


GetOrder(D::DerivativeOrder{O}) where {O} = O
GetOrder(O::Int) = O

Base.show(io::IO,DO::DerivativeOperator{TT,DIM,O}) where {TT,DIM,O} = print("Order ",GetOrder(DO.order)," ",DIM," dimensional diffusion SBP operator.")



# _next(i,start,stop) = i+1 > stop ? rem(i+1,stop+1)+start : i+1
# _prev(i,stop) = i-1 ≤ start ? mod1(i-1,stop) : i-1
@inline _next(i,stop) = i+1 > stop ? rem(i+1,stop+1)+2 : i+1
@inline _prev(i,stop) = i-1 ≤ 1 ? mod1(i-1,stop-1) : i-1

# Base.show(io::IO,DO::DerivativeOrder{O}) where O = print(io, "order ",O," second derivative operator")
