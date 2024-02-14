#  O seguinte algoritmo é dado pelo artigo intitulado "A Scaled Dai–Yuan Projection-Based Conjugate Gradient Method for Solving Monotone Equations with Applications"  


using LinearAlgebra

function merito(x,F)
    return  0.5* norm(Fx)^2
end

function CG_DYP(x, F; ϵ= 1.e-10, σ=1.e-4,β=0.6,m=0.8,n=-0.1, xu=0.8)
    k = 0
    xk = x          # xk means x_{k-1} in the article
    newx = x        # newx means x_k in the article
    # αk = b
    α= β
    Fx= F(x)
    newFx= Fx
    sk = newx-xk
    uk = sk
    # tk1= m *()
    while true
        if norm(F(xk)) > ϵ
            if k == 0
                dk = -F(xk) # Alteração, pois no artigo não mostra quem é d0 .  
            else 
                sk = newx-xk
                uk = sk
                yk = newFx - Fx
                pk = 2*(merito(xk,F)-merito(newx,F)) + sk'*(Fk+newFk) 
                wk = yk+ (xu * max(pk, 0)/ norm(sk)^2)*uk
                tk1 = m * (norm(wk)^2 / sk'*wk) - n * (sk'*wk / norm(sk)^2)
                dk = - Fk + ((wk-tk1*sk)'*Fk / dk'*wk)*dk
                if norm(dk) > ϵ
                    while -F(xk + α*dk)'*dk < σ* α* norm(dk)^2
                        α = α*β
                    end 
                    zk= xk+ αk*dk
                    if norm(zk) > ϵ
                        xk = newx
                        newx = newx - (F(zk)'*(newx - zk)/ norm(zk)^2)* F(zk)
                        k= k+1
                    else 
                        return k, xk, F(xk), norm(F(xk))     
                    end
                         
                else 
                    return k, xk, F(xk), norm(F(xk))   
                end
            end 
        else  
            return k, xk, F(xk), norm(F(xk)) 
        end
    end  
end

