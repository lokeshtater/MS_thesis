import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

data_2_pt = np.genfromtxt('./data_SR/tight_binding/2_pt_alt_100_sites.csv', delimiter=',', dtype=complex, skip_header=1) #each row contains data for one time-point
#data_4_pt = np.genfromtxt('./data_SR/alt_f.csv', delimiter=',', dtype=complex, skip_header=1) #each row contains data for one time-point


N= int(math.sqrt(len(data_2_pt[0])-1)) #Number of lattice sites


ind_D_mm = [] #column indices for elements D_mm  
# ind_F_mnmn = [] #column indices for elements F_mnmn
ind_D_mn = [] #column indices for elements D_ij

M = np.linspace(0,int(N/2 - 1),int(N/2)) #domain of fluctuation, half the system

for n in M:
    ind_d_mm  = int(N*n + n + 1)                  # "+1" because first column belongs to time
    ind_D_mm.append(ind_d_mm)
    for m in M:
        # ind_f_mnmn = int((N**3)*n + (N**2)*m + (N)*n + m + 1)
        # ind_F_mnmn.append(ind_f_mnmn)
        ind_D_mn.append(int(N*n +m +1))


#arrays to store time ordered data- \sum_{m,n \in M} D_mm(t), F_mnmn(t), |D_mn(t)|^2
time = []
D_mm = []
#F_mnmn = []
D_mn_sq = []


for row in data_2_pt:
    time.append(row[0])
    
    d_mm = 0.0 + 0.0j
    d_mn_sq = 0.0 + 0.0j
    
    for i in ind_D_mm: d_mm = d_mm + row[i]
    for i in ind_D_mn: d_mn_sq = d_mn_sq + abs(row[i])**2
    
    D_mm.append(d_mm)
    D_mn_sq.append(d_mn_sq)


    
# for row in data_4_pt:
#     f_mm = 0.0 + 0.0j
#     for i in diag_F: f_mm = f_mm + row[i]
#     F_mnmn.append(f_mm) 

T = range(len(time))

#w = [(-F_mnmn[t]+D_mm[t]-(D_mm[t]**2)) for t in T] #SR in terms of summations of 2 and 4 pt. correlators (general expression)

w_sq = [(D_mm[t]-D_mn_sq[t]) for t in T] #SR in terms of summations of 2_pt correlators, when Wick's thm is applicable


#saving SR into file:

filename = "data_SR/tight_binding/SR_alt_100_sites.csv"

df = pd.DataFrame(w_sq)
df.to_csv(filename,index=False)


