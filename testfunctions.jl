function probI(x)
    v = Float64[]
    push!(v,exp(x[1])-1)
    for i in 2:length(x)
        push!(v, exp(x[i])+x[i-1]-1)
    end
    return v
end

function projI(x)
    v = Float64[]
    for i in 1:length(x)
        push!(v, max(0,x[i]))
    end
    return v
end

