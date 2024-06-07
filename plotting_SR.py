import numpy as np
import math
import matplotlib.pyplot as plt

D1 = np.genfromtxt('./data_SR/alt_d.csv', delimiter=',', dtype=complex, skip_header=1) #each row contains data for one time-point
#F = np.genfromtxt('./data_SR/alt_f.csv', delimiter=',', dtype=complex, skip_header=1) #each row contains data for one time-point
D2 = np.genfromtxt('./data_SR/alt_f.csv', delimiter=',', dtype=complex, skip_header=1) #each row contains data for one time-point


N= int(math.sqrt(len(D1[0])-1)) #Number of lattice sites
diag_D = [] #coloumn for diagonal elements (D_mm) in the data  
diag_F = [] #coloumn for diagonal elements (F_mnmn) in the data
ind_ij = []
M = np.linspace(0,49,50) #domain of fluctuation

for i in M:
    diag_d = int(N*i + i + 1)                  # "+1" because first column belongs to time
    diag_D.append(diag_d)
    for j in M:
        diag_f = int((N**3)*i + (N**2)*j + (N)*i + j + 1)
        diag_F.append(diag_f)
        ind_ij.append(int(N*i +j +1))


time = []
D_mm = []
F_mnmn = []
D_ij_sq = []

D_mm_1 = []
D_ij_sq_1 = []


for row in D1:
    time.append(row[0])
    d_mm = 0.0 + 0.0j
    d_ij_sq = 0.0 + 0.0j
    for i in diag_D: d_mm = d_mm + row[i]
    for i in ind_ij: d_ij_sq = d_ij_sq + abs(row[i])**2
    
    D_mm.append(d_mm)
    D_ij_sq.append(d_ij_sq)

for row in D2:
    d_mm = 0.0 + 0.0j
    d_ij_sq = 0.0 + 0.0j
    for i in diag_D: d_mm = d_mm + row[i]
    for i in ind_ij: d_ij_sq = d_ij_sq + abs(row[i])**2
    
    D_mm_1.append(d_mm)
    D_ij_sq_1.append(d_ij_sq)

    
# for row in F:
#     f_mm = 0.0 + 0.0j
#     for i in diag_F: f_mm = f_mm + row[i]
#     F_mnmn.append(f_mm) 

T = range(len(time))

#w = [(-F_mnmn[t]+D_mm[t]-(D_mm[t]**2)) for t in T]
w = [(D_mm[t]-D_ij_sq[t]) for t in T]


#plt.plot(time,F_mnmn, marker='o', linestyle='-', color='b', label='F_mnmn')
#plt.plot(time,D_mm, marker='o', linestyle='-', color='g', label='D_mm')
# plt.plot(time,w, marker='o', linestyle='-', color='r', label='w')
# w1 = [(D_mm_1[t]-D_ij_sq_1[t]) for t in T]
# plt.plot(time,w1, marker='o', linestyle='-', color='b', label='w1')

# plt.xscale('log',base=10)


# plt.xlabel('time')
# plt.title('Plot of X vs Y')
# Add a legend
plt.legend()

plt.show()
