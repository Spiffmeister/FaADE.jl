using Test

using Pkg
Pkg.activate(".")


Dx = [-1.0,1.0]
Dy = [-π,π]
nx = 11
ny = 21



Dom = FaADE.Grid1D(Dx,nx)



Dom = FaADE.Grid2D(Dx,Dy,nx,ny)




function IdMap(X,x,p,t)
    X[1] = 0.0
    X[2] = 0.0
end

FaADE.construct_grid(IdMap,Dom,)


