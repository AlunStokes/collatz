using Printf
using Distributed

addprocs(4)

function collatz(x)
    l = Int64[]
    push!(l,x)
    while x > 1
        println(l)
        if x % 2 == 0
            x = x ÷ 2
        else
            x = 3*x + 1
        end
        push!(l, x)
    end
    return l
end

@everywhere D = Dict()

@everywhere function isCSelfCont(x)
    d = x
    while x > 1
        if haskey(D, x) == true
            D[d] = D[x]
            return D[x]
        end
        if x % 2 == 0
            x = x ÷ 2
        else
            x = 3*x + 1
        end
        if x % d == 0
            D[d] = true
            return true
        end
    end
    D[d] = false
    return false
end

function getCSelfCont(N)
    l = Int64[]
    n = 1
    while n <= N
        println(D)
        if isCSelfCont(n)
            push!(l, n)
        end
        n += 1
    end
    return l
end

#=
@time begin
    println(getCSelfCont(10^1))
end
=#

for n in UnitRange{UInt64}(1,10)
    println(n)
end