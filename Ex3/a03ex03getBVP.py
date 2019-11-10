import numpy as np
import scipy.sparse
import seaborn as sns
import matplotlib.pyplot as plt

def analytical_u(x):
	return 1+4*x**2-3*x**3

def analytical_f(x):
	return -3*x**3+40*x**2-14*x-7

def a03ex03getBVP(p):

	N=2**p-1
	h=1/(N+1)

	dd=scipy.sparse.diags([1,-2,1],[-1,0,1],shape=[N,N])
	d0=scipy.sparse.diags([-1/2,0,1/2],[-1,0,1],shape=[N,N])
	i=scipy.sparse.diags([0,1,0],[-1,0,1],shape=[N,N])

	l=-1/(h**2)*dd-4/h*d0+i

	x=1/h*np.asarray(range(1,N+1))
	f=analytical_f(x)

	f[0]=f[0]+(-2*h+1)/h**2
	f[-1]=f[-1]+(2*h+1)/h**2
	return x,l,f


def a03ex03solve():

	x,l,f=a03ex03getBVP(p)	
	u=analytical_u(x)
	u_h=scipy.sparse.linalg.spsolve(l, f) 
	error=max(abs(u-u_h))
	print(u_h)

	return error


error_list=[]
h_list=[]
for p in range(3,4):
	error=a03ex03solve()
	error_list.append(error)
	h_list.append(1/(2**p))
	print(error)

ax=sns.scatterplot(np.log(h_list),np.log(error))
ax.set(xlabel='log(h)', ylabel='log(error)')
plt.show()


