# O seguinte algoritmo é dado pelo artigo intitulado "Descent Perry conjugate gradient methods for systems of monotone nonlinear equations"


using LinearAlgebra, Plots
function CG_DY(x, F; ϵ, σ,β)
    k = 0
    xk = x   
    newx = x
    # αk = b
    Fx = F(x)
    fx = 1/2 * norm(Fx)^2


    if norm(F(xk)) > ϵ
        if k == 0
            dk = -F(xk) # Alteração, pois no artigo não mostra quem é d0 .  
        else
        yk = newFk - Fk
        sk = newx-xk
        uk = sk
        pk = 2*(fk-newfk) + sk'*(Fk-newFk)
        e=0.8
        wk = 
        m=0.8
        n=-0.1
        #tk1= m *(norm())
        end 
    else
        return k, xk, F(xk), norm(F(xk))
    end 

    
end