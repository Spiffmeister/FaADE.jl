
using Test
using LinearAlgebra
using FaADE




function buildprob2D(D,f,cx,cy,e)
    ue  = zeros(eltype(D.gridx),size(D))
    u   = zeros(eltype(D.gridx),size(D))
    kx  = zeros(eltype(D.gridx),size(D))
    ky  = zeros(eltype(D.gridx),size(D))
    Dxky = zeros(eltype(D.gridx),size(D))
    Dykx = zeros(eltype(D.gridx),size(D))
    for i in 1:D.nx
        for j in 1:D.ny
            if FaADE.Grid.coordtype(D) == CartesianMetric
                ue[i,j] = e(D.gridx[i],D.gridy[j])
                u[i,j]  = f(D.gridx[i],D.gridy[j])
            else
                ue[i,j] = e(D.gridx[i,j],D.gridy[i,j])
                u[i,j]  = f(D.gridx[i,j],D.gridy[i,j])
            end
        end
    end
    uxx = zeros(eltype(D.gridx),(D.nx,D.ny))
    return uxx, u, [kx, ky, Dxky, Dykx], ue
end

function setKoefficient!(K,kx,ky,G::FaADE.Grid.LocalGridType{TT,2,MET}) where {TT,MET}
    if MET == CartesianMetric
        for i in 1:G.nx
            for j in 1:G.ny
                K[1][i,j] = kx(G[i,j]...)                  # Kx
                K[2][i,j] = ky(G[i,j]...)                  # Ky
            end
        end
    else
        for i in 1:G.nx
            for j in 1:G.ny
                K[1][i,j] = kx(G[i,j]...) * G.J[i,j] * (G.qx[i,j]^2 + G.qy[i,j]^2)                  # Kx
                K[2][i,j] = ky(G[i,j]...) * G.J[i,j] * (G.rx[i,j]^2 + G.ry[i,j]^2)                  # Ky
                K[3][i,j] = kx(G[i,j]...) * G.J[i,j] * (G.qx[i,j]*G.rx[i,j] + G.qy[i,j]*G.ry[i,j])  # DxKy
                K[4][i,j] = ky(G[i,j]...) * G.J[i,j] * (G.qx[i,j]*G.rx[i,j] + G.qy[i,j]*G.ry[i,j])  # DyKx
            end
        end
    end
    K
end


# New solver 1 volume
cbottom(u) = [u,0.0]
cleft(v) = [0.0,v]
cright(v) = [1.0,v]
ctop(u) = [u,1.0]

Dom = Grid2D(cbottom,cleft,cright,ctop,21,21)

DomO = Grid2D([0.0,1.0],[0.0,1.0],21,21)

u₀(x,y) = x^2
# u₀(x,y) = exp.(-((x-0.5)^2 + (y-0.5)^2) / 0.02)
kx(x,y) = 2.0
ky(x,y) = 2.0

uxx, u, K, ue = buildprob2D(Dom,u₀,kx,ky,(x,y)->2.0)

setKoefficient!(K,kx,ky,DomO)


order = FaADE.Derivatives.DerivativeOrder{2}()


DOc = FaADE.Derivatives.DerivativeOperator{Float64,2,typeof(order),:Variable}(order,Dom.nx,Dom.ny,Dom.Δx,Dom.Δy,false,false)
DO = FaADE.Derivatives.DerivativeOperator{Float64,2,typeof(order),:Constant}(order,Dom.nx,Dom.ny,Dom.Δx,Dom.Δy,false,false)

FaADE.Derivatives.mul!(uxx,u,K,DOc)

FaADE.Derivatives.mul!(uxx,u,[K[1],K[2]],DO)



"""
    Test with packing x points near the origin
"""

u₀(x,y) = x^3

exact(x,y) = 6x

kx(x,y) = 1.0
ky(x,y) = 1.0

order = FaADE.Derivatives.DerivativeOrder{2}()

f(u) = (sinh(π*(2u-1))/sinh(π) + 1)/2 # point packing at 0.5


# New solver 1 volume
cbottom(u) = [f(u),0.0]
cleft(v) = [0.0,v]
cright(v) = [1.0,v]
ctop(u) = [f(u),1.0]

DomO = Grid2D([0.0,1.0],[0.0,1.0],11,11)
uxx, u, K, ue = buildprob2D(DomO,u₀,kx,ky,(x,y)->2.0)
setKoefficient!(K,kx,ky,DomO)

DO = FaADE.Derivatives.DerivativeOperator{Float64,2,typeof(order),:Constant}(order,DomO.nx,DomO.ny,DomO.Δx,DomO.Δy,false,false)
FaADE.Derivatives.mul!(uxx,u,[K[1],K[2]],DO)



DomC = Grid2D(cbottom,cleft,cright,ctop,11,11)
uxxc, u, K, ue = buildprob2D(DomC,u₀,kx,ky,(x,y)->2.0)
setKoefficient!(K,kx,ky,DomC)


DOc = FaADE.Derivatives.DerivativeOperator{Float64,2,typeof(order),:Variable}(order,DomC.nx,DomC.ny,DomC.Δx,DomC.Δy,false,false)
FaADE.Derivatives.mul!(uxxc,u,K,DOc)








