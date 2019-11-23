import numpy as np
import scipy.sparse as sp

def a05ex01_get_laplace(L,N,is_in_domain):
    # generate coordinates and grid spacing
    xx,yy = np.meshgrid(np.linspace(0,L,N), np.linspace(0,L,N))
    xx = xx.T
    yy = yy.T
    h  = L/(N-1)
    
    
    # determine points in bar_omega_h and lexicographical enumeration
    bar_omega_h   = is_in_domain(xx,yy,L)
    ix            = np.zeros((N,N))
    k             = np.sum(bar_omega_h.flatten(order="F"))
    ix[bar_omega_h] = np.arange(0,k)
    
    ## students part with some suggestions / hints
    #
    # Please follow the suggestions in the assignment to fill in the proper
    # functionality here. 
    
    omega_h = bar_omega_h # was meinst du mit return?
    gamma_h = bar_omega_h # 
    
    # sparse matrix assembly: indices ii,jj and values aa
    #
    # Q: What is ix_xmy,ix_xpy,... ?
    #
    #ii = np.concatenate((ix(omega_h);ix(omega_h);ix(omega_h);ix(omega_h);ix(omega_h);ix(gamma_h)))
    #jj = np.concatenate((ix(omega_h);ix_xmy(omega_h);ix_xpy(omega_h);ix_xym(omega_h);ix_xyp(omega_h);ix(gamma_h)))
    #
    # Q: what is k1,k2 ?
    #
    #aa = np.concatenate((np.ones(k1) *(+4)/h**2, # u(x,y)
    #      np.ones(k1) *(-1)/h**2, # u(x-1,y)
    #      np.ones(k1) *(-1)/h**2, # u(x+1,y)
    #      np.ones(k1) *(-1)/h**2, # u(x,y-1)
    #      np.ones(k1) *(-1)/h**2, # u(x,y+1)
    #      np.ones(k2)))
    #
    #Lh = sparse(ii(:),jj(:),aa(:));
    
    Lh = sp.csr_matrix((k,k))
    
    return  Lh,xx,yy,bar_omega_h,omega_h,gamma_h