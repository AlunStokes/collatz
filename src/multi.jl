using Distributed
n_cores = tryparse(Int, ARGS[2])
addprocs(n_cores)
@everywhere using SharedArrays

@everywhere function isCSelfCont(x)
    d = x
    while x > 1
        if x % 2 == 0
            x = x ÷ 2
        else
            x = 3*x + 1
        end
        if x % d == 0
            return true
        end
    end
    return false
end

@everywhere function getCSelfCont(L, U, offset, skip)
    n = Int64.(L) + Int64.(2 * offset + 1)
    if L % 2 == 1
        n = Int64.(L) + Int64.(2 * offset)
    end
    while n <= U
        if isCSelfCont(n)
            println(n)
        end
        n = n + 2 * skip
    end
end

@time begin
    k = tryparse(Int,ARGS[1])
    n_mach = tryparse(Int, ARGS[3])
    mach_index = tryparse(Int, ARGS[4])
    print("Checking up to 10^")
    println(k)
    #10^(10):10^(12) |> dump
    @sync @distributed for n in 0:(n_cores - 1)
    #println("I'm worker $(myid()), working on i=$i")
        #println(n)
        getCSelfCont(1, Int64.(10)^(k), n * n_mach + mach_index, n_cores * n_mach)
    end
end