from scipy import linalg
import scipy
from scipy.sparse import diags
import timeit
import numpy as np



def a01e03sparse(n):

	'''
	c) In cases where most entries of a matrix are zero, as is the case with tridiagonal matrices 
	and many other types, the sparse matrix format that assumes a value of zero for all 
	unspecified entries has a memory footprint that grows as O(n) and not O(n^2) as is the 
	case with dense matrices.
	'''
	return diags([-1,2,-1],[-1,0,1],shape=(n,n))


# b)



code1='scipy.sparse.linalg.spsolve(Kn,b)'

code2='np.linalg.solve(Kn.todense(),b)'

code3='np.linalg.inv(Kn.todense()).dot(b)'

codelist=[code1,code2,code3]


for n in (10,100,1000):


	setup='''
import scipy.sparse.linalg
import numpy as np
import random

def a01e03sparse(n):
	return scipy.sparse.diags([-1,2,-1],[-1,0,1],shape=(n,n))

n={}
random.seed(a=42)
b=np.random.rand(n,1)
Kn=a01e03sparse(n)
	'''.format(n)

	for i,code in enumerate(codelist):
		print('code {}'.format(i))
		print(timeit.timeit(setup=setup,stmt =code, number = 10))