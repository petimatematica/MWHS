# The presented algorithm is made using the article titled "Descent Perry conjugate gradient methods for systems of monotone nonlinear equations". Such an algorithm is done using the first given descent direction. 



# Function to execute the main algorithm
function algorithm1(F, x0)
    # Initialize variables
    iter = 0
    fn = 0
    x = x0
    F_km1 = NaN
    x_km1 = NaN
    d = NaN
    X, Y, Z = Float64[], Float64[], Float64[]
    

     t = time()
    # Infinite loop to iterate until convergence or max iterations
    while true
        Fx = F(x) # Compute F at x
        fn +=1
        norm_Fx = norm(Fx,2)  # Compute the 2-norm of F(x)


         # Print the current iteration and the norm of F(x)
        println("iter $iter  norm_Fx $norm_Fx")

        # Check for convergence
        if norm_Fx < epsilon
             et = time() - t
            return x,0,fn,et # Return current x and a success code
        end
       
        # Descent direction
        if iter == 0
            d = -Fx # Initial descent direction
        else
            y = Fx - F_km1 # Difference of F(x) between iterations
            s = x - x_km1 # Step size
            f_k = merit1(x,F)    # Compute merit function at current x
            f_km1 = merit1(x_km1,F) # Compute merit function at previous x
            zeta = 2.0*(f_km1 - f_k) + dot(s,F_km1 + Fx)  # Compute zeta
            mu = s
            w = y + xi * (max(zeta,0.0) / dot(s,mu)) * mu  # Compute w
            stw = dot(s,w)
            t1 = (1 + (m_star * norm(w,2)^2 / stw - n_star * stw / norm(s,2)^2)) #Compute t2
            dtw = dot(d,w)
            d = -Fx + (dot(w-t1*s,Fx) / dtw) * d    # Update descent direction
        end

        # Linesearch to find the step size alpha
        alpha, iter_int = linesearch1(x,d,F)
        fn += iter_int

        # Compute new candidate solution
        z = x + alpha * d

        F_z = F(z) #Compute F at new canditate solution
        fn += 1
        norm_F_z = norm(F_z,2) #Compute the 2-norm of F(z) 

        # Check for convergence 
        if norm_F_z < epsilon
             et = time()-t
            return z,2,fn,et # Return new solution and a success code
        end

        # Update variables for the next iteration
        x_km1 = x
        
        x = x - (dot(F_z,x-z) / (norm_F_z)^2) * F_z  #Projection step 

        F_km1 = Fx

        iter += 1 # Increment iteration counter

        # Check for maximum iterations
        if iter > maxiter
             et = time()- t
            return x,1,fn,et # Return current x and a failure code
        end
        
    end

end

# Function to perform Linesearch
function linesearch1(x,d,F)
    
    m = 0
    fn_linesearch=0
    snorm_d2 = sigma * norm(d,2)^2 # Precompute sigma * ||d||^2

    # Infinite loop to find the suitable step size alpha
    while true
        alpha = beta^m #Compute alpha
        if alpha < 1.e-6
            return  1.e-6,fn_linesearch   # Return minimum alpha if below threshold
        end
      #  println("alpha = $alpha")
        q = x + alpha * d # Compute new candidate solution
        F_q = F(q) # Compute F at the new candidate solution
        fn_linesearch +=1

        # Check the stopping test condition
        stptest = dot(F_q,d) + alpha * snorm_d2 > 0.0
        if ~stptest
            return alpha,fn_linesearch # Return alpha if stopping test is satisfied
        end
        m += 1 # Increment the counter
    end


end

function merit1(x,F)

    return 0.5 * norm(F(x))^2

end