using LinearAlgebra, Plots
function CG_DY(x, F; ϵ, σ,β)
    k = 0
    xk = x
    newx = x
    # αk = b
    f = F(x)


    if norm(F(xk)) > ϵ
        if k == 0
            dk = -F(xk) # Alteração, pois no artigo não mostra quem é d0 .  
        else
            
        end 
    else
        return k, xk, F(xk), norm(F(xk))
    end 

    
end