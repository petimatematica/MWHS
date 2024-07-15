#
#
using LinearAlgebra

include("MCGP1.jl")
#include("MCGP2.jl")
include("testfunctions.jl")
# Dimension seting
# n = 50

# initial guess
function PontosIniciais(n)

    x1 = Float64[]
    x2 = Float64[]
    x3 = Float64[]
    x4 = Float64[]
    x5 = Float64[]
    x6 = Float64[]
    x7 = Float64[]
    x8 = Float64[]

    for i in 1:n
        push!(x1,1)
        push!(x2,2-(i/n))
        push!(x3,1/i)
        push!(x4,(i/n)-1)
        push!(x5,(-1)^i*0.25)
        push!(x6,i/n)
        push!(x7,i^(-1))
        push!(x8,1-i/n)
    end

    return x1, x2, x3, x4, x5, x6, x7, x8
end

x1, x2, x3, x4, x5, x6, x7, x8 = PontosIniciais(1000)

# Code parameters

maxiter = 1000
epsilon = 1.e-10
sigma = 1.e-4
beta = 0.6
xi = 0.8
m_star = 0.8
n_star = -0.1
q_star = 0.8
r_star = -0.5
sol = algorithm1(prob1, x2);
 #sol = algorithm1(prob8, x2);

 println("$(sol[3])   $(sol[4])")


