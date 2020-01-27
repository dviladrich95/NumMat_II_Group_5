import numpy as np
import scipy.sparse as sps
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
np.set_printoptions(precision=2)

def f_function(x):
	return 0*np.sin(2*np.pi*x)+0

def g1_function(x):
	return 0*np.sin(2*np.pi*x)+0

def g2_function(x):
	return 2*np.cos(x)


#N1=N_r N2=N_phi
def a05ex02solvePDE(N_r,N_phi):

	h_phi=2*np.pi/(N_phi+1)
	phi_i=0
	phi_f=2*np.pi-h_phi

	h_r=1/(N_r+1)
	r_i=1+h_r
	r_f=2-h_r

	phi=np.linspace(phi_i,phi_f,num=N_phi+1)
	r=np.linspace(r_i,r_f,num=N_r)

	i_coef=1/(h_r**2)*sps.diags([1,-2,1],[-1,0,1],shape=[N_r,N_r])

	r_m1=sps.diags(1/(r[1:N_r]),-1,shape=[N_r,N_r])
	r_p1=sps.diags(1/(r[0:N_r-1]),1,shape=[N_r,N_r])
	r_coef=1/(2*h_r)*(r_p1-r_m1)

	rr_coef=sps.diags(1/(r*r),0,shape=[N_r,N_r])

	t=sps.diags([1,-2,1],[-1,0,1],shape=[N_phi+1,N_phi+1]).tolil()
	t[0,N_phi]=1
	t[N_phi,0]=1
	t=1/(h_phi**2)*t

	dd_r=sps.kron(i_coef,sps.eye(N_phi+1))
	d_r=sps.kron(r_coef,sps.eye(N_phi+1))
	dd_phi=sps.kron(rr_coef,t)

	l=-(dd_r+d_r+dd_phi)



	g1_r=np.zeros(N_r)
	g1_r[0]=1
	g1=np.kron(g1_r,g1_function(np.ones(N_phi+1)))

	g2_r=np.zeros(N_r)
	g2_r[N_r-1]=1
	g2=np.kron(g2_r,g2_function(phi))

	f=1/(h_r**2)*(g1+g2)
	f=f.transpose()

	#print(g2_function(phi))
	#print(g2)
	#print(f)
	u_h=sps.linalg.spsolve(l, f)
	print(l.shape)

	u_h=np.reshape(u_h,(N_r,N_phi+1))
	return u_h
	#rr,phiphi=np.meshgrid(r,phi)

	#t is the laplace matrix


def exact_solution(r,phi):

	angle=2*np.cos(phi)	
	radial=2/3*(r-1/r)
	u=angle*radial
	return u


for q in range(3,15):
	N_r=2**q
	N_phi=2**q-1

	print(N_r,N_phi)

	h_phi=1/(N_phi+1)
	phi_i=h_phi
	phi_f=2*np.pi-h_phi

	h_r=1/(N_r+1)
	r_i=1+h_r
	r_f=2-h_r

	phi=np.linspace(phi_i,phi_f,num=N_phi+1)
	r=np.linspace(r_i,r_f,num=N_r)
	PHIPHI,RR=np.meshgrid(phi,r)
	u=exact_solution(RR,PHIPHI)

	u_h=a05ex02solvePDE(N_r,N_phi)

	ax=sns.heatmap(u_h,xticklabels=False,yticklabels=False)
	ax.set(xlabel='phi', ylabel='r')

	plt.savefig('a06ex02c_uh.png')
	plt.clf()
	ax=sns.heatmap(u-u_h,xticklabels=False,yticklabels=False)
	ax.set(xlabel='phi', ylabel='r')
	plt.savefig('a06ex02c_error.png')
	max_error=np.max(u-u_h)
	print(max_error,'max error')

#ax.set(xlabel='x', ylabel='u_h')


'''
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
'''