
####
####
# Objective function section   "of the article "Descent Perry conjugate gradient methodsfor systems of monotone nonlinear equations" 


#Problem 1  
function prob1(x)
    v = Float64[]
    n = length(x)
    for i in 1:n
        push!(v, log(x[i]+1)- (x[i]/n))   
    end
    return v
end


#Problem 2 Já ok
function prob2(x)
    v = Float64[]
    n = length(x)
    for i in 1:n
        push!(v,2*x[i] - sin(abs(x[i])))
    end
    return v
end


#Problem 3  VERIFICAR

function prob3(x)
    v = Float64[]
    n = length(x)
    push!(v, 2*x[1]+ sin(x[1]) - 1)
  
    for i in 2:n-1
       push!(v, -2*x[i-1]+ 2*x[i]+ sin(x[i]) - 1)
    end
    push!(v, 2*x[n] + sin(x[n])-1)
    return v
end

 
#Problem 4  VERIFICARRR

function prob4(x)
    v = Float64[]
    n = length(x)
    c = 0.999
    mu(j) = (j-0.5)/n
    for i in 1:n
        push!(v , x[i]- (1 - (c/(2*n))* sum(mu(i) * x[j]/(mu(i)+mu(j)) for j in 1:n))^-1)
    end
    return v
end

#Problem 5 
function prob5(x)

    v = Float64[]
    n = length(x)

    push!(v, 2.5*x[1] + x[2] - 1)

    for i in 2:n-1
       push!(v, x[i-1]+ 2.5*x[i]+ x[i+1] - 1)
    end

     push!(v, x[n-1] + 2.5*x[n] - 1)
    return v
end

#Problem 6 
function tridgen1(a,b,c,n)
    di=a*ones(n-1)                  # di - diagonal inferior
    d=b*ones(n)                     # d - diagonal
    ds=c*ones(n-1)                  # ds - diagonal superior
    return Tridiagonal(di, d, ds)
end


function prob6(x) 
    n = length(x)
    b1 = Float64[]

    for i in 1:n
        push!(b1,exp(x[i])-1)
    end

    A1 = tridgen1(-1,2,-1,n)

    return A1*x+b1
end


#Problem 7 

function tridgen2(a,b,c,n)
    di=a*ones(n-1)                  # di - diagonal inferior
    d=b*ones(n)                     # d - diagonal
    ds=c*ones(n-1)                  # ds - diagonal superior
    return Tridiagonal(di, d, ds)
end


function prob7(x)
    n = length(x)
    b2 = Float64[]

    for i in 1:n
        push!(b2,sin(x[i])-1)
    end

    A2 = tridgen2(0,2,-1,n)

    return A2*x+b2
end


#Problem 8  

function tridgen3(a,b,c,n)
    di=a*ones(n-1)                  # di - diagonal inferior
    d=b*ones(n)                     # d - diagonal
    ds=c*ones(n-1)                  # ds - diagonal superior
    return Tridiagonal(di, d, ds)
end


function prob8(x)
    n = length(x)
    b3 = Float64[]

    for i in 1:n
        push!(b3,-n)
    end

    A3 = tridgen3(2,5,3,n)

    return A3*x+b3
end


#Problem 9 

function prob9(x)
    v = Float64[]
    n = length(x)

    push!(v, -2*x[1]- x[2] + exp(x[1]) - 1)

    for i in 2: n-1
       push!(v, -x[i-1] + 2*x[i] - x[i+1] + exp(x[i]) - 1)
    end 
        push!(v, 2*x[n]- x[n-1] + exp(x[n])- 1)  
    
    return v
    
end

#Problem 10 

function prob10(x)
    v = Float64[]
    n = length(x)
    for i in 1:n
        push!(v, x[i] - 1/n*x[i]^2 + 1/n * sum( x[i] + i for i in i:n) )
    end
    return v
end

#Problem 11 

function prob11(x)
    v = Float64[]
    n = length(x)
    push!(v, x[1] - exp(cos((x[1]+x[2])/(n+1))))

    for i in 2:n-1
        push!(v, x[i] - exp(cos((x[i-1] + x[i] + x[i+1])/(n+1))) )
    end

    push!(v, x[n] - exp(cos((x[n-1] + x[n])/(n+1))))  
  
     return v
end



#Problem 12 

function prob12(x)
    v = Float64[]
    n = length(x)
    push!(v, x[1]*(x[1]^2 + x[2]^2) - 1 )

    for i in 2: n-1
        push!(v, x[i]*(x[i-1]^2 + 2*x[i]^2 + 2*x[i+1]^2)- 1)
     end
         push!(v, x[n]*(x[n-1]+x[n]^2))  
    
     return v
end

