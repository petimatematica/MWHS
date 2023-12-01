using LinearAlgebra, Plots

function CG_DY(x, F, proj; a = 0.1, g = 0.0001, c = 0.99, b = 1, ϵ = 1.e-9, λ_min = 1, λ_max = 1)
    k = 0
    xk = x
    newx = x
    αk = b
    if norm(F(xk)) > ϵ # Alteração de: norm(F(xk)) > c
        d0 = -F(xk)  
    else
        return k, xk, F(xk), norm(F(xk))
    end

    if dot(-F(xk + αk*d0),d0) ≥ g*αk*norm(F(xk + αk*d0))*norm(d0)^2
        if c > 1
            while dot(-F(xk + αk*d0),d0) ≥ g*αk*norm(F(xk + αk*d0))*norm(d0)^2
                αk = αk*c 
            end
            αk = αk/c
        end
    end

    if norm(F(xk + αk*d0)) == 0
        return k, xk, F(xk), norm(F(xk))
    else
        vk = dot(F(xk + αk*d0),-αk*d0)/norm(F(xk + αk*d0))^2
        newx = proj(xk-vk*F(xk + αk*d0))
    end

    k = k + 1

    while true
        if norm(F(xk)) > ϵ # Alteração de: norm(F(xk)) > c
            hk = xk + αk*d0
            sk = hk-xk
            yk = F(hk) - F(xk) + a*sk
            σk = dot(yk-sk,F(newx))/norm(newx)^2
            σk = min(1,abs(σk))
            τk = 1+σk*dot(F(newx),d0)/dot(d0,yk)
            βk = norm(F(newx))^2/dot(d0,yk) # Alteração do algoritmo proposto
            d0 = -τk*F(newx) + σk*βk*d0
            
            αk = b

            if dot(-F(newx + αk*d0),d0) ≥ g*αk*norm(F(newx + αk*d0))*norm(d0)^2
                if c > 1
                    while dot(-F(newx + αk*d0),d0) ≥ g*αk*norm(F(newx + αk*d0))*norm(d0)^2
                        αk = αk*c 
                    end
                    αk = αk/c
                end
            end

            if norm(F(newx + αk*d0)) == 0
                return k, newx, F(newx), norm(F(newx))
            else
                vk = dot(F(newx + αk*d0),-αk*d0)/norm(F(newx + αk*d0))^2
                xk = newx
                newx = proj(newx-vk*F(newx + αk*d0))
            end
        
            k = k + 1
        else
            return k, newx, F(newx), norm(F(newx))
        end
    end
end