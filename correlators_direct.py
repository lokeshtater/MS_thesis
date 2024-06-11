import numpy as np #basic packages
import time
from scipy.linalg import expm
import cmath
import pandas as pd 

N = 125 #number of lattice sites, sites labelled from 0 to N-1


#initial states for 2-point and 4-point correlator vectors:

D0 = np.zeros(N**2,dtype=complex) #vector to store initial state
#F0 = np.zeros(N**4,dtype=complex)

for n in range(N):
    for m in range(N):
        if(n==m):
           if((n+1)%2==0):
               D0[N*n + m] = 1.0 + 0.0j     #alternating pure state 
          #if(n<int(N/2)):             #domain wall state 
            #    D0[N*n + m] = 1.0 + 0.0j      
            # if(n<int(N/2)):                  #staggered state
            #     if((n+1)%2==0):
            #         D0[N*n + m] = 1.0 + 0.0j
        
print("Initial State Defined")

#defining evolution matrices for 2-point and 4-point correlator vectors:
A = np.zeros((N,N),dtype=complex)
Id = np.zeros((N,N),dtype=complex)

#tight-binding
# for i in range(N):
#     for j in range(N):
#         if(i==j):  Id[i][j] = 1.0 #also defining identity matrix
#         if (abs(i-j)==1):  A[i][j] = -1.0j
# A[0][N-1] = A[N-1][0] = -1.0j #periodicity of lattice

#long range hopping: #closed boundaries
alpha = 1.4

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
df1.to_csv(f'./data_SR/long_hopping/2_pt_alt_{N}_sites.csv')


print(df1.head())
print("DONE!")

