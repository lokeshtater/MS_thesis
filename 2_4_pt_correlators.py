import numpy as np
import pandas as pd
import time
from joblib import Parallel, delayed

from qutip import *


import Operators
from scipy.integrate import solve_ivp


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
D_evol = np.kron(Id,A) + np.kron(B,Id)
F_evol = np.kron(Id,np.kron(Id,np.kron(Id,A))) +  np.kron(Id,np.kron(Id,np.kron(A,Id))) +  np.kron(Id,np.kron(B,np.kron(Id,Id))) +  np.kron(B,np.kron(Id,np.kron(Id,Id)))

print("Evolution Matrices Defined")

####Solving for evolution of 2 and 4 point correlator vectors- D and F#############
###################################################################################

#Select which ODE to solve for
y0 = D0   
evol_matrix = D_evol


def ode_func(t,y): #function to define ODE
    return np.dot(evol_matrix,y)
    
    
t_span = (0,10**5)



print("Solving ODE...")
#solve ODE
start_time = time.time()
solution = solve_ivp(ode_func, t_span, y0, t_eval=t_eval, args=(evol_matrix,))
end_time = time.time()
print(f"ODE solved. Time taken:{end_time-start_time}")



t_points = solution.t
y_points = solution.y.T  # Transpose to get a 2D array with time points along rows

# Create a DataFrame to store the results
df = pd.DataFrame(data=y_points, index=t_points)
df.index.name = 'time'

# Save the DataFrame to a CSV file
df.to_csv('./SR_data/alt_d.csv')


print("DONE!")
print(df.head())  # Print the first few rows to verify



