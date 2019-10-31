import math
import numpy as np

def a02e02solveheat(mode,x,t,a,c):

	#dirichlet b. c.
	if mode==1:
		n=np.asarray(range(1,100,2))
		sine=np.sin(math.pi*n*x/math.sqrt(a))	

		b=16.0/(math.pi*(4*n-n**3)) 
		
		g=b*sine
		f=np.exp(-((math.pi*n)**2+c)*t)

		u=np.sum(g*f)

	#von neumann b. c.
	elif mode == 2:
		cosine=1-np.cos(2*math.pi*x/math.sqrt(a))	

		f=np.exp(-((math.pi*2)**2+c)*t)

		u=np.sum(cosine*f)


	return u

x=np.arange(0,1,0.1)
t=np.arange(0,1,0.1)

xx,tt=np.meshgrid(x, t)

np.meshgrid
a02e02solveheat(1,xx,tt,1,1)



