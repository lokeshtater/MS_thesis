import numpy as np
from qutip import *
import Operators

N = 2 #number of lattice sites, sites labelled from 0 to N-1
a = Operators.annihilation(N) #list of single site annihilation operators for lattice of size N
vacc =  tensor([basis(2,0) for i in range(N)]) #vaccum state


#defining initial staggered state
stagg = vacc
for i in range(N):
    if i%2==0:
        stagg = a[i].dag() * stagg

init = stagg

D = np.empty((N, N), dtype=complex)
F = np.empty((N, N, N, N), dtype=complex)

for m in range(N):
    for n in range(N):
        D[m,n] = expect(a[m].dag()*a[n],init)
        for p in range(N):
            for q in range(N):
                F[m,n,p,q] = expect(a[m].dag()*a[n].dag()*a[p]*a[q],init)
        
        


