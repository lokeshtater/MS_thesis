import numpy as np
import math
import matplotlib.pyplot as plt

data_2_pt = np.genfromtxt('./data_SR/tight_binding/2_pt_alt_75_sites.csv', delimiter=',', dtype=complex, skip_header=1) #each row contains data for one time-point
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
F_mnmn = []
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

w_sq = [(D_mm[t]-D_mn_sq[t]) for t in T] #SR in terms of summations of 2_pt correlators, when Wick's thm holds true



#fitting

#first linear region: #region chosen visually and is t<<1/J
###############################
t1 = 10**(-3)
t2 = 10**(-1) 

time_1 = [t for t in time if t1<t<t2]

log_t = np.log10(time_1)
indices = [i for i,t in enumerate(time) if t1<t<t2]
log_w_sq = np.log10([w_sq[i] for i in indices])

beta1 = np.polyfit(log_t,log_w_sq,1)[0] #fitting log(w^2) against log(t)

w_sq_fit1 = [t**beta1 for t in time_1]

#second linear region: #region chosen visually and is t~1/J
###############################
t1 = 10**(0)
t2 = 10**(1) 

time_2 = [t for t in time if t1<t<t2]

log_t = np.log10(time_2)
indices = [i for i,t in enumerate(time) if t1<t<t2]
log_w_sq = np.log10([w_sq[i] for i in indices])

beta2 = np.polyfit(log_t,log_w_sq,1)[0] #fitting log(w^2) against log(t)

w_sq_fit2 = [t**beta2 for t in time_2]


##Plotting####

plt.plot(time,w_sq, marker='o', linestyle='-', color='r', label=f'N={N}') #raw data
plt.plot(time_1,w_sq_fit1, linestyle='-', color='b', label=fr"$t^{{{round(beta1,3)}}}$") #fit of region1 (t<<1/J)
plt.plot(time_2,w_sq_fit2, linestyle='-', color='g', label=fr"$t^{{{round(beta2,3)}}}$") #fit of region2 (t~1/J)



plt.xscale('log',base=10)
plt.yscale('log',base=10)

plt.xlabel('Jt (J = nearest neighbour hopping stregth)')
plt.ylabel(r"Bipartite particle number fluctuations")
plt.title('Free fermions, tight binding, periodic boundary conditions')

plt.legend()
plt.show()
