import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from a02e02solveheat import a02e02solveheat

def a02e02plot():

	x=np.arange(0,1,0.01)
	t=np.arange(0,10,0.01)
	log_t=np.asarray([0,0.001,0.01,0.1,1])

	xx,tt=np.meshgrid(x, log_t)

	#xx,tt=np.meshgrid(x, t)
	#sns.heatmap(a02e02solveheat(2,xx,tt,1,1))

	for mode in [1,2]:

		for c in [0,1]:

			with sns.cubehelix_palette(8, start=2, rot=0, dark=0, light=.95, reverse=True):
				for heat in a02e02solveheat(mode,xx,tt,1,c):
					sns.scatterplot(x,heat)
				plt.savefig('a1_mode{}_c{}'.format(mode,c))
				plt.show()

	for mode in [1,2]:
		try:
			with sns.cubehelix_palette(8, start=2, rot=0, dark=0, light=.95, reverse=True):
				for heat in a02e02solveheat(mode,xx,tt,-1,0):
					sns.scatterplot(x,heat)
				plt.savefig('a-1_mode{}_c0'.format(mode))
				plt.show()
		except:
			print('Negative thermal diffusivity a')

'''
When trying to calculate the function with a negative thermal coefficient a the equation fails
when trying to compute the square root of a, which is set to be a strictly positive parameter
at the beginning because changing the sign of the parameter a, known in physics as the thermal 
diffusivity, also changes the fundamental behaviour of the PDE, making it an ill posed initial
value problem.
'''


a02e02plot()