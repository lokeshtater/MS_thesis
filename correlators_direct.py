import numpy as np
import pandas as pd
import time
from scipy.linalg import expm
import cmath

from qutip import *


import Operators


N = 12 #number of lattice sites, sites labelled from 0 to N-1
a = Operators.annihilation(N) #list of single site annihilation operators for lattice of size N
vacc =  tensor([basis(2,0) for i in range(N)]) #vaccum state


#defining pure initial staggered state
stagg = vacc
for i in range(N):
    if (i+1)%2==0:  #all odd sites filled
        stagg = a[i].dag() * stagg

init = stagg

D0 = np.empty(N**2,dtype=complex)
F0 = np.empty(N**4,dtype=complex)

#initial states for 2-point and 4-point correlator vectors:
for n in range(N):
    for m in range(N):
        D0[N*n + m] = expect(a[m].dag()*a[n],init) 
        for p in range(N):
            for q in range(N):
                F0[q*N**3 + p*N**2 + n*N + m] = expect(a[m].dag()*a[n].dag()*a[p]*a[q],init)

print("Initial State Defined")

#defining evolution matrices for 2-point and 4-point correlator vectors:
A = np.zeros((N,N),dtype=complex)
Id = np.zeros((N,N),dtype=complex)

for i in range(N):
    for j in range(N):
        if(i==j):  Id[i][j] = 1.0 #also defining identity matrix
        if (i==((j+1)%N) | i==((j-1)%N) ):  A[i][j] = -1.0j
       

B = -A

#Evolution
T = 10**6 #time upto which we evolve
times = np.linspace(0,int(T),100)


D = np.empty((len(times),int(N**2)),dtype=complex)
F = np.empty((len(times),int(N**4)),dtype=complex)

for i in range(len(times)):
    start = time.time()
    t = times[i]
    D[i] = np.dot(np.kron(expm(B*t),expm(A*t)),D0)
    F[i] = np.dot(np.kron(expm(B*t),np.kron(expm(B*t),np.kron(expm(A*t),expm(A*t)))),F0)
    
    end = time.time()
    print(f"iteration:{i} time:{t}  time taken:{end-start} seconds")

    
df1 = pd.DataFrame(data=D, index=times)
df1.index.name = 'time'

df2 = pd.DataFrame(data=F, index=times)
df2.index.name = 'time'

# Save the DataFrame to a CSV file
df1.to_csv('./data_SR/alt_d.csv')
df2.to_csv('./data_SR/alt_f.csv')

print("DONE!")
#print(df.head())  # Print the first few rows to verify

