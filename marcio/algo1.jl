function algorithm1(F)
    iter = 0
    x = copy(x0)
    while true
        Fx = F(x)
        norm_Fx = norm(Fx,2)

        if norm_Fx < epsilon
            return x,0
        end
       
        # Descent direction
        if iter == 0
            d = -Fx
        else
            y = Fx - F_km1
            s = x - x_km1
            f_k = merit(x)
            f_km1 = merit(x_km1)
            zeta = 2.0*(f_km1 - f_k)+dot(s,F_km1 + Fx)
            mu = s
            w = y + xi * max(zeta,0.0) / dot(d,mu) * mu
            t = Parei aqui
            d = -Fx + ( dot())
        end

        # Linesearch
        alpha = linesearch(x,d,F)

        z = x + alpha * d

        F_z = F(z)
        norm_F_z = norm(F_z,2)

        if norm_F_z < epsilon
            return z
        end

        x_km1 = x
        
        x = x - dot(F_z,x-z) / norm_F_z^2 * F_z

        iter += 1
        if iter > maxiter
            return x,1
        end
        
        F_km1 = Fx
    end

end

# Linesearch
function linesearch(x,d,F)

    m = 0
    snorm_d2 = sigma * norm(d,2)^2
    while true
        alpha = beta^m
        q = x + alpha * d
        F_q = F(q)
        stptest = dot(F_q,d) + alpha * snorm_d2 > 0.0
        if ~stptest
            return alpha
        end
        m += 1
    end


end

function merit(x)

    return 0.5 * norm(F(x))^2

end