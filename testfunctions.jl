
####
####
# Objective function 

#Problem 1
function prob1(x)
    v = Float64[]
    n = length(x)
    for i in 1:n
        push!(v, log(x[i]+1)- (x[i]/n))   
    end
    return v
end


#Problem 2
function prob2(x)
    v = Float64[]
    for i in 1:length(x)
        push!(v,2*x[i] - sin(abs(x[i])))
    end
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



























# # Problems to numerical results


# #Problem 1

# function prob1(x)
#     v = Float64[]
#     for i in 2:length(x)
#         push!(v, log10(x[i]+1)- (x[i]/n))   
#     end
#     return v
# end

# #Problem 2

# function prob2(x)
#     v = Float64[]
#     for i in 1:length(x)
#         push!(v, 2*x[i]- sin(abs(x[i])))   
#     end
#     return v 
# end

# #Problem 3




# # function projI(x)
# #     v = Float64[]
# #     for i in 1:length(x)
# #         push!(v, max(0,x[i]))
# #     end
# #     return v
# # end

# #Primeiro problema para teste dado por : 
# # Problem I - F_1 (x)= e^{x_1} - 1, 
# #             F_i(x) = e^{x_i} + x{i-1} -1, i= 2,3,...,n-1, and  Ω = R_+^n
 

# # function prob1(x)
# #     v = Float64[]
# #     push!(v,exp(x[1])-1)
# #     for i in 2:length(x)
# #         push!(v, exp(x[i])+x[i-1]-1)
# #     end
# #     return v
# # end