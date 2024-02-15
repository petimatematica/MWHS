#
# 
#
using LinearAlgebra

include("algo1.jl")
# Dimension seting
n = 3

# Objective function
function prob2(x)
    v = Float64[]
    for i in 1 : n
        push!(v,2*x[i] - sin(abs(x[i])))
    end
    return v
end

# initial guess
x0 = ones(n)

# Code parameters
maxiter = 1000
epsilon = 1.e-8
sigma = 0.75
beta = 0.5
xi = 1.5
sol = algorithm1(prob2)


println(sol)