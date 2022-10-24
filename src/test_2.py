import numpy as np
import math
import cmath
import time
from numpy import linalg as LA
import functools

def testFun( a, b, c ):
	return a*b+c
	
sigy = np.array( [[0,-1j],[1j,0]] )
vec1 = np.array( [1,1j] )

vFun = np.vectorize( testFun )

ln = 100000
a = np.random.uniform (high = 1, size = ln )
b = np.random.uniform( high = 1, size = ln )
c = np.random.uniform( high = 1, size = ln )

result1 = np.zeros( ln )
tic = time.time()
for ita in range(ln):
			result1[ita] = testFun( a[ita], b[ita], c[ita] )
toc = time.time()
time1 = toc - tic

tic = time.time()
result2 = vFun( a, b, c )
toc = time.time()
time2 = toc - tic

tic = time.time()
result3 = testFun( a, b, c )
toc = time.time()
time3 = toc - tic

print( time1, time2, time3 )

	



print( vFun( [1,2,3], [4,5,6], [7,8,9] ) )

# U_mat = np.array( [[-0.02,0.95],[0.95,+0.02]] )

# eValLst, eVecLst = LA.eig(U_mat)
# eVecLstTr = np.transpose(eVecLst)

# print(U_mat)
# print(eValLst)
# print(eVecLst)
# print(eVecLstTr)

# file_arr=open('arr_data.csv','ab')
# np.savetxt(file_arr, [eValLst], delimiter=',',fmt='%f')

# for n in range(0,N,1):
	# axs[n].plot( vLst[1], chernLst[indCenter,:,n] )
	# axs[n].set( ylabel= ("C"+str(n)) )
	# chernSum += chernLst[indCenter,:,n]

N = 5
ang_divide = 20
ang2pi = 2 * math.pi
ang_ln = ang_divide + 1
ang_step = ang2pi / ang_divide
ang_lst = ang_step * np.array( list( range(ang_ln) ) )
print( ang_lst )