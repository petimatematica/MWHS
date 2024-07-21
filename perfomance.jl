

using LinearAlgebra, DataFrames, Plots, BenchmarkProfiles

include("MCGP1.jl")
include("MCGP2.jl")
include("testfunctions.jl")

#Parâmetros

maxiter = 10000
epsilon = 1.e-6
sigma = 1.e-4
beta = 0.6
xi = 0.8
m_star = 0.8
n_star = -0.1
q_star = 0.8
r_star = -0.5

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
        push!(x2, n^2)
        push!(x3,1/i)
        push!(x4,1-(i/n))
        push!(x5,(-1)^i*0.25)
        push!(x6,i/n)
        push!(x7,i/(n-1))
        push!(x8,1-i/n)
    end

    return x1, x2, x3, x4, x5, x6, x7, x8
end

pontos_iniciais = Float64[]

dim= [500]
estrategias= ["MCGP1", "MCGP2"]
func=[prob2, prob3, prob4, prob5, prob5, prob6, prob7, prob8, prob9, prob10, prob11, prob12]
#, prob4, prob5, prob5, prob6, prob7, prob8, prob9, prob10, prob11, prob12
times1 = Float64[]
times2 = Float64[]
iter1 = Float64[]
iter2 = Float64[]
aval1 = Float64[]
aval2 = Float64[]

t_inicial = time()

for F in func
    for dimension in dim
        x1, x2, x3, x4, x5, x6, x7, x8 = PontosIniciais(dimension)
        pontos_iniciais = [x1, x2, x3, x4, x5, x6, x7, x8] 
        for x0 in pontos_iniciais
            # println(x0)
            
            for estrategia in estrategias
                println("Rodando teste com: F=$F, dim=$dimension, estrategia=$estrategia")
                t0= time()
                if estrategia == "MCGP1"
                    (x,ierror,fn,et,seqx) = algorithm1(F, x0)
                    println("ierror = ", ierror)
                    t_decorrido= time() - t0 

                    if ierror == 1
                        push!(times1, Inf)
                        push!(iter1, Inf)
                        push!(aval1, Inf)
                        else
                        push!(times1, t_decorrido)
                        push!(iter1, size(seqx, 2))
                        push!(aval1, fn)  
                    end
                    
                end
                if estrategia == "MCGP2"
                    (x,ierror,fn,et, seqx) = algorithm2(F, x0)
                    t_decorrido= time() - t0 
                    println("ierror = ", ierror)
                    if ierror == 1
                        push!(times2, Inf)
                        push!(iter2, Inf)
                        push!(aval2, Inf)
                        else
                        push!(times2, t_decorrido)
                        push!(iter2, size(seqx, 2))
                        push!(aval2, fn)  
                    end
                end
            end 
        end 
    end
end
t_final= time() - t_inicial
problems= length(pontos_iniciais)* length(dim)* length(pontos_iniciais)* length(func)
println("Número de problemas=", problems) 
println("Tempo total gasto =", t_final/60)

X= [times1 times2]
Y= [iter1 iter2]
Z= [aval1 aval2]

colors= [:blue2, :green1]


P1 = performance_profile(PlotsBackend(), X, ["MCGP1", "MCGP2"], 
xlabel = "CPU time ratio", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 2, dpi = 1000)

P2 = performance_profile(PlotsBackend(), Y, ["MCGP1", "MCGP2"], 
xlabel = "Iteration", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 2, dpi = 1000)

P3 = performance_profile(PlotsBackend(), Z, ["MCGP1", "MCGP2"], 
xlabel = "Function evaluations", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 2, dpi = 1000)

println("Performance profiles generated, saving figures...")

savefig(P1, "Time_teste_500.png")
savefig(P2, "Aval_500.png")
savefig(P3, "iter_500.png")

println("Figures saved, the code has finished running.")