import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def a02e01solution(r,phi):


	angle=2*np.cos(phi)	
	radial=2/3*(r-1/r)
	u=angle*radial

	return u

def a02e01plot():

	r=np.arange(1,2,0.01)
	phi=np.arange(0,2*math.pi,0.01)

	R,PHI=np.meshgrid(r, phi)

	x=np.arange(-2,2,0.005)
	y=np.arange(-2,2,0.005)
	
	X,Y=np.meshgrid(x, y)
	R=np.sqrt(X**2+Y**2)
	PHI=np.arctan2(Y,X)

	cmap=sns.palettes.diverging_palette(240,10,n=51)
	heat=a02e01solution(R,PHI)
	heat[R<1]=None
	heat[R>2]=None
	sns.heatmap(heat,cmap=cmap)
	plt.show()


a02e01plot()




