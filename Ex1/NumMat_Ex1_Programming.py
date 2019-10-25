from scipy import linalg
import scipy
from scipy.sparse import diags
import timeit
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg
import random

def a01e03sparse(n):

	'''
	c) In cases where most entries of a matrix are zero, as is the case with tridiagonal matrices 
	and many other types, the sparse matrix format that assumes a value of zero for all 
	unspecified entries has a memory footprint that grows as O(n) and not O(n^2) as is the 
	case with dense matrices. Additionally there are matrix solving algorithms that take
	advantage of the sparcity of the matrix to improve solution speed.
	Note that from the results of the time_ratio_matrix the better performance of the 
	sparse matrix solvers might only overtake general methos when dimensionality n is large.

	'''
	return diags([-1,2,-1],[-1,0,1],shape=(n,n)).tocsr()


# b)


code1='scipy.sparse.linalg.spsolve(Kn,b)'
code2='np.linalg.solve(Kn.todense(),b)'
code3='np.linalg.inv(Kn.todense()).dot(b)'

codelist=[code1,code2,code3]

time_ratio_matrix=np.zeros((3,3))
mean_residual_matrix=np.zeros((3,3))

for i,n in enumerate((10,100,1000)):


	setup='''
import scipy.sparse.linalg
import numpy as np
import random

def a01e03sparse(n):
	return scipy.sparse.diags([-1,2,-1],[-1,0,1],shape=(n,n)).tocsr()

n={}
random.seed(a=42)
b=np.random.rand(n,1)
Kn=a01e03sparse(n)
	'''.format(n)

	for j,code in enumerate(codelist):
		print('code {}'.format(j+1))
		time=timeit.timeit(setup=setup,stmt =code, number = 10)

		if j==0:
			first=time
		time_ratio_matrix[i][j]=time/first
		print(time)





for i,n in enumerate((10,100,1000)):

	random.seed(a=42)
	b=np.random.rand(n,1)
	Kn=a01e03sparse(n)


	x1=scipy.sparse.linalg.spsolve(Kn,b)
	x2=code2=np.linalg.solve(Kn.todense(),b)
	x3=code3=np.linalg.inv(Kn.todense()).dot(b)

	mean_residual_matrix[i][0]=linalg.norm(Kn.dot(x1)-b)
	mean_residual_matrix[i][1]=linalg.norm(Kn.dot(x2)-b)
	mean_residual_matrix[i][2]=linalg.norm(Kn.dot(x3)-b)


print('mean_residual_matrix\n',mean_residual_matrix)
print('time_ratio_matrix\n',time_ratio_matrix)