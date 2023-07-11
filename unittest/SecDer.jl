

using FaADE
using BenchmarkTools
# using ProfileView
# using Cthulhu
using Profile

function buildgrid(n)
    x = collect(range(0.0,stop=1.0,length=n+1))
    Δx = x[2]-x[1]
    return n+1, x, Δx
end


nx,x,Δx = buildgrid(100)
ny,y,Δy = buildgrid(100)

order = 4
# 1D
u = x.^2
uxx = zeros(nx)
k = ones(nx)
D₂!(uxx,u,k,nx,Δx,order,0.0)
println(uxx)

# @benchmark D₂!($uxx,$u,$k,$nx,$Δx,$order,0.0)

# D₂!(uxx,u,k,nx,Δx,order,0.0)
# Profile.clear_malloc_data()
# D₂!(uxx,u,k,nx,Δx,order,0.0)


# D1uxx = zeros(nx);
# D1u = x.^3;
# D1k = ones(nx);
# u = 6x;
# @benchmark D₂!($D1uxx,$D1u,$D1k,$nx,$Δx,2)
# D₂!(D1uxx,D1u,D1k,nx,Δx,2)
# function f(uxx,u,k,n,Δ,o)
#     for i = 1:1000000
#         D₂!(uxx,u,k,n,Δ,o,0.0)
#     end
# end

# @profview f(uxx,u,k,nx,Δx,2)
# @profview f(uxx,u,k,nx,Δx,2)




#=
# 2D
u = zeros(Float64,(nx,ny))
for i = 1:nx; for j = 1:ny
    u[i,j] = x[i]^4 #+ y[j]^4
end; end

k = ones(Float64,(nx,ny))
uxx = zeros(Float64,(nx,ny))
D₂!(uxx,u,k,k,nx,ny,Δx,Δx,2,2)
=#


# uxx1 = zeros(nx,ny)
# uxx2 = zeros(nx,ny)

# Profile.clear_malloc_data()
# D₂!(uxx,u,k,k,nx,ny,Δx,Δx,2,2)
# uxx'




# using BenchmarkTools
# @benchmark D₂!($uxx,$u,$k,$k,$nx,$ny,$Δx,$Δy,2,2)

# D₂!(uₓₓ::AbstractMatrix{T},u::AbstractMatrix{T},c::AbstractMatrix{T},n::Integer,Δ::T,dim::Integer=1,order::Integer=2)

# function DeeXX(uxx,uxx1,uxx2,u,k,nx,ny,Δx,Δy,o)
#     D₂!(uxx1,u,k,nx,Δx,1,o)
#     D₂!(uxx2,u,k,ny,Δy,2,o)
#     @. uxx = uxx1 + uxx2
# end


# @benchmark DeeXX($uxx,$uxx1,$uxx2,$u,$k,$nx,$ny,$Δx,$Δy,2)

