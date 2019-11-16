import numpy as np
import scipy.sparse as sps
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

def analytical_u(x):
	return 1+4*x**2-3*x**3

def analytical_f(x):
	return -3*x**3+40*x**2-14*x-7

def get_x_list(N):
	x=[]
	while len(x)!=N:
		x=np.random.uniform(0,N,size=N-2)
		print(x)
		x=np.insert(x,0,0)
		print(x)
		x=np.append(x,1)
		print(x)
		x=np.unique(x)
		print(x)
	return x

def a04ex02getPDE(x,f,consts,flag):
	
	a=consts[0]
	b=consts[1]
	c=consts[2]
	alpha=consts[3]
	beta=consts[4]
	h=np.diff(x)
	h=[x[i]-x[i-1] for i in range(1,len(x))]

	f_dd_m1= lambda i: 2/(h[i]*(h[i]+h[i+1]))
	f_dd_0= lambda i: -2/(h[i]*h[i+1])
	f_dd_p1= lambda i: 2/(h[i+1]*(h[i]+h[i+1]))

	f_d_0_1= lambda i: -1/(h[i]+h[i+1])
	f_d_0_2= lambda i: 1/(h[i]+h[i+1])

	f_d_p_1= lambda i: -1/h[i+1]
	f_d_p_2= lambda i: 1/h[i+1]

	f_d_m_1= lambda i: -1/h[i]
	f_d_m_2= lambda i: 1/h[i]


	N=len(h)-1

	dd_m1=[f_dd_m1(i) for i in range(1,N)]
	dd_0=[f_dd_0(i) for i in range(N)]
	dd_p1=[f_dd_p1(i) for i in range(N-1)]
	dd=np.diag(dd_m1,k=-1)+np.diag(dd_0,k=0)+np.diag(dd_p1,k=1)
	i=sps.eye(N)

	if flag=='-':
		d_1=[f_d_m_1(i) for i in range(1,N)]
		d_2=[f_d_m_2(i) for i in range(N)]

		d=np.diag(d_1,k=-1)+np.diag(d_2,k=0)
	elif flag=='+':
		d_1=[f_d_p_1(i) for i in range(N)]
		d_2=[f_d_p_2(i) for i in range(N-1)]
		d=np.diag(d_1,k=0)+np.diag(d_2,k=1)
	elif flag=='0':
		d_1=[f_d_0_1(i) for i in range(1,N)]
		d_2=[f_d_0_2(i) for i in range(N-1)]
		d=np.diag(d_1,k=-1)+np.diag(d_2,k=1)
	l=-a*dd+b*d+c*i

	if flag=='-':
		f[0]=f[0]+(a*f_dd_m1(0)-b*f_d_m_1(0))*alpha
		f[-1]=f[-1]+a*f_dd_p1(N-1)*beta
	elif flag=='+':
		f[0]=f[0]+a*f_dd_m1(0)
		f[-1]=f[-1]+(a*f_dd_p1(N-1)-b*f_d_p_2(N-1))*beta
	elif flag=='0':
		f[0]=f[0]+(a*f_dd_m1(0)-b*f_d_0_1(0))*alpha
		f[-1]=f[-1]+(a*f_dd_p1(N-1)-b*f_d_0_2(N-1))*beta
	
	return l,f

def a04ex02solve(x,f,consts,flag):

	l_h,f_h=a04ex02getPDE(x,f,consts,flag)
	u_h=sps.linalg.spsolve(l_h,f_h)
	
	alpha=consts[3]
	beta=consts[4]

	u_h=np.insert(u_h,0,alpha)
	u_h=np.append(u_h,beta)
	return u_h

'''
N=10
consts=[1,-4,1,1,2]
x=np.linspace(0,1,num=N+2)
ax=sns.scatterplot(x,analytical_u(x))
for noise in np.linspace(0,1,num=2):
	x=np.linspace(0,1,num=N+2)+noise/2*np.random.uniform(-1/(N+2),1/(N+2),N+2)
	x[0]=0
	x[-1]=1
	f=analytical_f(x)[1:-1]
	flag='0'
	u_h=a04ex02solve(x,f,consts,flag)
	ax=sns.scatterplot(x,u_h)
ax.set(xlabel='x', ylabel='u_h')
plt.show()
'''