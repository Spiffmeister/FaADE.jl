
using Test
using LinearAlgebra
using FaADE




function buildprob2D(D,f,cx,cy,e)
    ue  = zeros(eltype(D.gridx),size(D))
    u   = zeros(eltype(D.gridx),size(D))
    kx  = zeros(eltype(D.gridx),size(D))
    ky  = zeros(eltype(D.gridx),size(D))
    kxy = zeros(eltype(D.gridx),size(D))
    for i in eachindex(D.nx)
        for j in eachindex(D.ny)
            ue[i,j] = e(D.gridx[i],D.gridy[j])
            u[i,j]  = f(D.gridx[i],D.gridy[j])
        end
    end
    uxx = zeros(eltype(D.gridx),(D.nx,D.ny))
    return uxx, u, [kx, ky, kxy], ue
end

function setKoefficient!(K,kx,ky,G::FaADE.Grid.LocalGridType{TT,2,MET}) where {TT,MET}
    for i in 1:G.nx
        for j in 1:G.ny
            K[1][i,j] = kx(G[i,j]...) * G.J[i,j] * (G.qx[i,j]^2 + G.qy[i,j]^2)
            K[2][i,j] = ky(G[i,j]...) * G.J[i,j] * (G.rx[i,j]^2 + G.ry[i,j]^2)
            K[3][i,j] = kx(G[i,j]...) * G.J[i,j] * (G.qx[i,j]*G.rx[i,j] + G.qy[i,j]*G.ry[i,j])
        end
    end
    K
end


# New solver 1 volume
cbottom(u) = [u,0.0]
cleft(v) = [0.0,v]
cright(v) = [1.0,v]
ctop(u) = [u,1.0]

Dom1V = Grid2D(cbottom,cleft,cright,ctop,21,21)

u₀(x,y) = x^2
# u₀(x,y) = exp.(-((x-0.5)^2 + (y-0.5)^2) / 0.02)
kx(x,y) = 2.0
ky(x,y) = 2.0

uxx, u, K, ue = buildprob2D(Dom1V,u₀,kx,ky,(x,y)->2.0)

setKoefficient!(K,kx,ky,Dom1V)




