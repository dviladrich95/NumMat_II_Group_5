import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def a02e02solveheat(mode,x,t,a,c):

	#dirichlet b. c.
	if mode==1:
		def sine(n):
			sine=np.sin(math.pi*n*x/math.sqrt(a))
			return sine

		def b(n):
			b=16.0/(math.pi*(4*n-n**3)) 
			return b
		
		def g(n):
			return b(n)*sine(n)

		def f(n):
			return np.exp(-((math.pi*n)**2+c)*t)

		def u(n):
			return sum([g(i)*f(i) for i in range(1,n,2)])

		#n=100 sets a truncation of the fourier series that is computed fast while having 
		#a truncation error for (x,t)=(0.5,0) of around e-6, more than necessary for plotting
		#and with enough precision to use as input of other transforms.
		return u(100)

	#von neumann b. c.
	elif mode == 2:
		constant=1
		cosine=np.cos(2*math.pi*x/math.sqrt(a))	
		decay_constant=np.exp(-c*t)
		decay_cosine=np.exp(-((math.pi*2)**2+c)*t)
		u=constant*decay_constant-cosine*decay_cosine

	return u