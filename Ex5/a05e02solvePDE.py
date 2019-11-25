import numpy as np
import scipy.sparse as sps
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
np.set_printoptions(precision=3)

def f_zero(x):
	return 0*np.sin(2*np.pi*x)+0

def f_function(x):
	return 1*np.sin(2*np.pi*x)+1




#data=[c, g0, g1, alpha0,alpha1, beta0, beta1]
def a05ex02solvePDE(N,f,data,flag):

	c, g0, g1, alpha0,alpha1, beta0, beta1=data

	if flag==1:

		x=np.linspace(0,1,num=N+2)
		h=(x[N+1]-x[0])/(N+1)
		dd=1/(h**2)*sps.diags([1,-2,1],[-1,0,1],shape=[N+2,N+2])
		i=sps.eye(N+2)
		l=-dd+c*i

		b0=np.zeros(N+2)
		b1=np.zeros(N+2)

		if beta0==0:
			b0[0]=alpha0-c*h**2/2
			l[0]=b0
			f[0]=g0

		else:
			b0[0],b0[1]=(1-(h/beta0*alpha0-c*h**2/2))/h**2,-1/h**2
			l[0]=b0
			f[0]=f[0]/2+g0/(h*beta0)


		if beta1==0:
			b1[N+1]=alpha1-c*h**2/2
			l[N+1]=b1
			f[N+1]=g1

		else:
			b1[N],b1[N+1]=-1/h**2,(1-(h/beta1*alpha1-c*h**2/2))/h**2
			l[N+1]=b1
			f[N+1]=f[N+1]/2+g1/(h*beta1)

		if alpha0==0 and alpha1==0 and c==0:
			h_add=sps.csr_matrix.transpose(sps.csr_matrix(np.ones(N+2)))
			v_add=sps.csr_matrix(np.transpose(np.ones(N+3)))
			v_add[0,N+2]=0
			l=sps.hstack((l,h_add))
			l=sps.vstack((l,v_add))
			f=np.concatenate((f,[0]),axis=0)
			u_h,q=np.split(sps.linalg.spsolve(l, f),[N+2])
			print(q,'q')

		else:
			u_h=sps.linalg.spsolve(l, f)



	elif flag==2:

		x=np.linspace(0,1,num=N+1)
		h=(x[N]-x[0])/(N+1)
		dd=1/(h**2)*sps.diags([1,-2,1],[-1,0,1],shape=[N+1,N+1]).tocsr()

		dd[0,N]=1/h**2
		dd[N,0]=1/h**2

		i=sps.eye(N+1)
		l=-dd+c*i

		h_add=sps.csr_matrix.transpose(sps.csr_matrix(np.ones(N+1)))
		v_add=sps.csr_matrix(np.transpose(np.ones(N+2)))
		v_add[0,N+1]=0
		l=sps.hstack((l,h_add))
		l=sps.vstack((l,v_add))
		f[N+1]=0
		u_h,q=np.split(sps.linalg.spsolve(l, f),[N+1])
	return x,u_h




#c)
N=500
h=1/(N+1)
x=np.linspace(0,1,num=N+2)
f=f_zero(x)
#data=[c, g0, g1, alpha0,alpha1, beta0, beta1]
data=[0,-0.5,2,0,0,1,1]
flag=1
x,u_h_c=a05ex02solvePDE(N,f,data,flag)
ax=sns.lineplot(x,u_h_c)
ax.set(xlabel='x', ylabel='u_h')
plt.show()

#d)
N=2**10
h=1/(N+1)
x=np.linspace(0,1,num=N+2)
f=f_function(x)
#data=[c, g0, g1, alpha0,alpha1, beta0, beta1]
data=[0,0,0,0,0,1,1]
flag=1
x,u_h_final=a05ex02solvePDE(N,f,data,flag)

error_list=[]
h_list=[]
for p in range(3,10):
	N=2**p
	h=1/(N+1)
	h_list.append(h)
	x=np.linspace(0,1,num=N+2)
	f=f_function(x)
	#data=[c, g0, g1, alpha0,alpha1, beta0, beta1]
	data=[1,0,0,0,0,0,0]
	flag=2
	x,u_h=a05ex02solvePDE(N,f,data,flag)
	error_list.append(sum(abs(u_h_final[::2**(10-p)]-u_h)))
	with sns.cubehelix_palette(8, start=2, rot=0, dark=0, light=.95, reverse=True):
		ax=sns.lineplot(x,u_h)

error_list=np.log(np.asarray(error_list))
h_list=np.log(np.asarray(h_list))
slope, intercept, r_value, p_value, std_err=stats.linregress(h_list,error_list)

ax.set(xlabel='x', ylabel='u_h')
plt.show()
print(slope,'slope')
ax=sns.scatterplot(h_list,error_list)
ax.set(xlabel='log(h)', ylabel='log(error)')
plt.show()