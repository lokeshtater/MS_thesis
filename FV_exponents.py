import numpy as np
import matplotlib.pyplot as plt

lat_size = [75,100] #lattice sizes, taken from original code
time = [(10**i) for i in np.linspace(-3,+3,100)] #time points of data, taken from original code

#dictionaries
w_sq = {}  #for raw SR data
beta1 = {} #growth1 exponent
w_sq_fit1 = {} 
beta2 = {}
w_sq_fit2 = {}
w_sq_sat_val = {} #averaged saturation value

t1,t2 = 10**(-3),10**(-1) #region1 (growth1)
t3,t4 = 10**(0),10**(1) #region2  (growth2)
t5,t6 = 10**(2),10**(3) #region3 (saturation)

indices1 = [i for i,t in enumerate(time) if t1<t<t2] 
indices2 = [i for i,t in enumerate(time) if t3<t<t4] 
indices3 = [i for i,t in enumerate(time) if t5<t<t6] 
time1 = [time[i] for i in indices1]
time2 = [time[i] for i in indices2]
time3 = [time[i] for i in indices3]

colors_blue = {50:"#d4f1f4",75:"#75e6da",100:"#189ab4",125:"#05445e"}
colors_green = {50:"#61d095",75:"#48bf84",100:"#439775",125:"#2a4747"}
colors_orange = {50:"#f5bb00",75:"#ec9f05",100:"#d76a03",125:"#bf3100"}

for N in lat_size:
    w_sq[N] = np.genfromtxt(f"./data_SR/tight_binding/SR_dom_{N}_sites.csv", delimiter=',', dtype=complex, skip_header=1)

    #fitting in region1

    log_w_sq = np.log10([w_sq[N][i] for i in indices1])
    log_t = np.log10(time1)
    beta1[N] = np.polyfit(log_t,log_w_sq,1)[0] #fitting log(w^2) against log(t)
    
    w_sq_fit1[N] = [t**beta1[N] for t in time1]
    

    #fitting in region2

    log_w_sq = np.log10([w_sq[N][i] for i in indices2])
    log_t = np.log10(time2)
    beta2[N] = np.polyfit(log_t,log_w_sq,1)[0] #fitting log(w^2) against log(t)
    
    w_sq_fit2[N] = [t**beta2[N] for t in time2]
    
    
    plt.plot(time,w_sq[N],marker='o', color=colors_blue[N], label=f'N={N}') #raw data
    plt.plot(time1,w_sq_fit1[N], linestyle='-', color=colors_green[N], label=fr"$t^{{{round(beta1[N],3)}}}$") #fit of region1 (t<<1/J)
    plt.plot(time2,w_sq_fit2[N], linestyle='-', color=colors_orange[N], label=fr"$t^{{{round(beta2[N],3)}}}$") #fit of region2 (t~1/J)



plt.xscale('log',base=10)
plt.yscale('log',base=10)

plt.xlabel('Jt (J = nearest neighbour hopping stregth)')
plt.ylabel(r"$w^{{2}}_{L/2}(L,t)$")
plt.title('Free fermions, tight binding, periodic boundary conditions')

plt.legend()
plt.show()
