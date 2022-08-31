





# struct SAT end
# function ()(uₓₓ,)









function boundary(order,Δ)
    
    BDxT = BDₓᵀ(order,Δ)

    α,τ = SATpenalties()
    if forcing 
        α = -α; τ = -τ
    else
    end
    
    SAT_Dirichlet_internal!(SAT,node,u,c,Δ,α,τ,BDxT,order)
end









