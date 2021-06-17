using Distributed
n_cores = tryparse(Int,ARGS[3])
addprocs(n_cores)
@everywhere using SharedArrays

@everywhere function isNegInv(x, k)
    r = x - invmod(k, x)
    d = x
    while x > 1
        if x % 2 == 0
            x = x ÷ 2
        else
            x = 3*x + 1
        end
        if x % d == r
            return true
        end
    end
    return false
end

@everywhere function getCSelfCont(U, k, n_cores, offset)
    n = Int64.(k + 1 + offset)
    while n <= U
        if n % k == 0
            n = n + n_cores
            continue
        end
        if isNegInv(n, k)
            println(n)
        end
        n = n + n_cores
    end
end

@time begin
    k = tryparse(Int,ARGS[1])
    a = tryparse(Int, ARGS[2])
    print("Checking up to 10^")
    print(a)
    print(" for k = ")
    println(k)
    #10^(10):10^(12) |> dump
    @sync @distributed for n in 0:(n_cores - 1)
    #println("I'm worker $(myid()), working on i=$i")
        #println(n)
        getCSelfCont(Int64.(10)^(a), k, n_cores, n)
    end
end