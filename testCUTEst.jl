

### 

using CUTEst, NLPModels, LinearAlgebra, DataFrames, Random, Printf
include("MCGP1.jl")
include("MCGP2.jl")
# Defina uma matriz de nomes de problemas CUTEst
problems = ["PALMER5D"]
dimension = ["2"]

for i in 1:length(problems)

    # Inicialize um vetor de modelos CUTEst
    nlp = CUTEstModel(problems[i])

    # Initial guess from CUTEst
    x0 = nlp.meta.x0

    # Objective functions

    global function f(x)
        return obj(nlp,x)  # Objective function from CUTEst
    end

    # Gradient of Objective function

    global function gradf(x)
        return grad(nlp,x) # Gradient of objective function (CUTEst)
    end



  # Code parameters
  global maxiter = 100000
  global epsilon = 1.e-5
  global sigma = 1.e-4
  global beta = 0.6
  global xi = 0.8
  global m_star = 0.8
  global n_star = -0.1
  global q_star = 0.8
  global r_star = -0.5
  global sol = algorithm1(gradf,x0)
  
  #println("$sol")

finalize(nlp)
end