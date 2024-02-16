#
# 
#
using LinearAlgebra

include("algo1.jl")
# Dimension seting
n = 2000

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
xi = 0.8
m_star = 0.8
n_star = -0.1
sol = algorithm1(prob2)


println(sol)