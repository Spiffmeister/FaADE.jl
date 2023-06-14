using Test

using Pkg
Pkg.activate(".")


Dx = [-1.0,1.0]
Dy = [-π,π]
nx = 11
ny = 21



Dom = SBP_operators.Grid1D(Dx,nx)



Dom = SBP_operators.Grid2D(Dx,Dy,nx,ny)




function IdMap(X,x,p,t)
    X[1] = 0.0
    X[2] = 0.0
end

SBP_operators.construct_grid(IdMap,Dom,)


