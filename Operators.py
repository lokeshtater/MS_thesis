from qutip import *

#defining destruction operator for i-th site
def single_destruct(m,N):
    site_ops = []
    for i in range(N):
        if (i==m):
            site_ops.append(destroy(2))
        else:
            site_ops.append(qeye(2))

    return tensor(site_ops)

#returns a site ordered list of annihilation operators: 
def annihilation(N):
    return [single_destruct(i,N) for i in range(N)]

