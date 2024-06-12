import numpy as np #basic packages
import time
from scipy.linalg import expm
import cmath
import pandas as pd 
import sys

N = int(sys.argv[1]) #number of lattice sites, sites labelled from 0 to N-1
system = sys.argv[2]
state = sys.argv[3]

print(f"{N} {state}")
#initial states for 2-point and 4-point correlator vectors:

D0 = np.zeros(N**2,dtype=complex) #vector to store initial state
#F0 = np.zeros(N**4,dtype=complex)


if(state=="alt"): ind = [(N*n + n) for n in range(N) if (n+1)%2==0]
if(state=="dom"): ind = [(N*n + n) for n in range(int(N/2))]
if(state=="stagg"): ind = [(N*n + n) for n in range(int(N/2)) if (n+1)%2==0]

D0[ind] = 1.0 + 0.0j

print("Initial State Defined")

#defining evolution matrices for 2-point and 4-point correlator vectors:
A = np.zeros((N,N),dtype=complex)
Id = np.zeros((N,N),dtype=complex)


if(system=="tight_binding"):
    for i in range(N):
        Id[i][i] = 1.0
        A[i][int((i+1)%N)] = A[i][i-1] = -1.0j

if(system=="long_hopping"):
    alpha = float(sys.argv[4])
    for i in range(N):
        Id[i][i] = 1.0
        for j in range(i+1,N):
            A[i][j] = A[j][i] = (1/(j-i)**alpha)*1.0j

B = -A


#Evolution######
#################

T = 3 #evolution will be computed for times 10^-T to 10^T
times = [(10**i) for i in np.linspace(-T,+T,100)]


D = np.zeros((len(times),int(N**2)),dtype=complex)
#F = np.zeros((len(times),int(N**4)),dtype=complex)


for i in range(len(times)):
    start = time.time()
    
    t = times[i]
    D[i] = np.dot(np.kron(expm(B*t),expm(A*t)),D0)
    #F[i] = np.dot(np.kron(expm(B*t),np.kron(expm(B*t),np.kron(expm(A*t),expm(A*t)))),F[i-1])
        
        
    end = time.time()
    print(f"iteration:{i} time:{t}  time taken:{end-start} seconds")

    
df1 = pd.DataFrame(data=D, index=times)
df1.index.name = 'time'
df1.to_csv(f'./data_SR/{system}/2_pt_{state}_{N}_sites.csv')


print(df1.head())
print(alpha)
print("DONE!")

