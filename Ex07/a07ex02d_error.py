import numpy as np
import scipy.sparse as sps
import seaborn as sns
import matplotlib.pyplot as plt

from a07ex02d import a07ex02getsold
from a07ex02c import a07ex02c


N=10
Nx=N
Ny=N
Nz=N
R=10
M=10

error_list=[]
t_list=np.linspace(0.1,0.5,30)

for T in t_list:

    ud=a07ex02getsold(T,Nx,Ny,Nz,R)
    uc=a07ex02c(N,M,T)
    udiff=ud-uc
    max_diff=np.max(np.abs(udiff))
    error_list.append(max_diff)
    print(max_diff)

error_list=np.asarray(error_list)
sns.scatterplot(t_list,np.log(error_list-0.0007092680976878638))
plt.show()