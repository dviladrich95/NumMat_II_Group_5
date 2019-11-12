import numpy as np
import scipy.sparse
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

def analytical_u(x):
	return 1+4*x**2-3*x**3

def analytical_f(x):
	return -3*x**3+40*x**2-14*x-7

def get_x_list(N):
	x=[]
	while len(x)!=N
		x_=np.random.uniform(size=N-2)
		x.insert(0,0)
		x.insert(-1,1)
		x=numpy.unique(x)
	return x

def get_h_list(x):


def a04ex02getPDE(x):
	h=numpy.diff(x)
	N=len(h)-1
	dd_m1=[2/(h[i]*(h[i]+h[i+1])) for i in range(N)]
	dd_0=[-2/(h[i]*h[i+1]) for i in range(N)]
	dd_p1=[2/(h[i+1]*(h[i]+h[i+1])) for i in range(N)]

	d0_m1=[-1/(h[i]+h[i+1]) for i in range(N)]
	d0_p1=[1/(h[i]+h[i+1]) for i in range(N)]

	dd=a*(np.diags(dd_m1)+np.diags(dd_0)+np.diags(dd_p1))
	d0=b*(np.diags(d0_m1)+np.diags(d0_p1))
	i=c*(scipy.sparse.eye(N)

	l=-1/(h**2)*dd-4/h*d0+i
	x=h*np.asarray(range(1,N+1))
	f=analytical_f(x)
	f[0]=f[0]+(-2*h+1)/h**2
	f[-1]=f[-1]+(2*h+2)/h**2
	return x,l,f


def a03ex03solve():

	x,l,f=a03ex03getBVP(p)	
	u=analytical_u(x)
	u_h=scipy.sparse.linalg.spsolve(l, f) 
	error=max(abs(u-u_h))
	return error


error_list=[]
h_list=[]
for p in range(2,15):
	error=a03ex03solve()
	error_list.append(error)
	h_list.append(1/(2**p))


slope, intercept, r_value, p_value, std_err=stats.linregress(np.log(h_list),np.log(error_list))
print('convergence rate:',slope)

ax=sns.scatterplot(np.log(h_list),np.log(error_list))
ax.set(xlabel='log(h)', ylabel='log(error)')
plt.show()




