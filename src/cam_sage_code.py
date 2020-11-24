from sage.all import *

#Functions sans documentation (sorry!)
def collatz(n):
    oplist=[]
    while n > 1:
        if n%2==0:
            n=n/2
            oplist.append(0)
        else:
            n=(3*n+1)/2
            oplist.append(1)
    oplist.append(1)
    return(oplist)

def revcollatz(L):
    n=1
    M=len(L)
    for jj in range(M-1):
        if L[M-2-jj] == 0:
            n=2*n
        else:
            n= (2*n-1)/3
    return(n)

# The following snippet tests that revcollatz is left inverse to collatz
# for the first bunch of integers
for n in range(1,100000):
    if revcollatz(collatz(n)) != n:
        print(n,"Error!")

# this snippet computes the sequence of collatz operations, stopping the first time you hit 1
# adjust k to get more data
k=10
data=[]
for n in range(1,2^k+1):
    LL=collatz(n)
    data.append([n,LL,len(LL)])
show(table(data))


# this snippet computes sequnce of revcollatz ops on the binary expansions of integers
# up to 2^k
k=10
data=[]
for n in range(1,2^k+1):
    L=ZZ(n).digits(2)
    data.append([L,revcollatz(L)])
show(table(data))
