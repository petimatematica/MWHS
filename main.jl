include("search.jl"); include("testfunctions.jl")

iter, x, Fx, normFx = CG_DY(ones(500), probI, projI, ϵ= 1.e-9) 
