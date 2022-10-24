import numpy as np
import math
import cmath
from numpy import linalg as LA
import functools

def PairLst(list1, list2): 
      
    paired_list = tuple(zip(list1, list2))  
    return paired_list 

def SortPairLst( pairLst ):
	def getKey(item):
		return item[0]
	return sorted(pLst,key=getKey)
	
def cmplxZeros( shape ):
	return np.zeros(shape,dtype=np.complex_)

def testFun( a, b ):
	return a+b

b=3

testFunA = functools.partial( testFun, b=b )

print("test numbers")
n = 2
print(n)
mat0 = np.zeros( (3,3) )
mat1 = np.ones( (3,3) )
matId = np.identity(3)
mat2 = matId;
mat2[1][1] = 2;
matProduct = matId.dot(mat1)
print("mat2",mat2);
print("Product",matProduct)

im_num = cmath.exp(1j * 2 * math.pi)
print(im_num)

print("cos", math.cos( math.pi / 3 ))
print(np.ones(5))

for ind in range(1,6):
	print(ind)

print(4 * matId)

ang2pi = 2 * math.pi
ang_divide = 10
ang_ln = ang_divide + 1
ang_step = ang2pi / ang_divide
ang_lst = np.zeros(ang_ln);
for ind in range(ang_ln):
	ang_lst[ind] = ind * ang_step
print(ang_lst)

w, v = LA.eig(np.diag((1, 2, 3)))
print(w)
print(v)

print(cmath.log(1j))

print(np.log([1,3,-2+1j]))

#pLst = PairLst([2,1,3],[['b','b'],['a','a'],['c','c']])
pLst = PairLst([2,1,3],['b','a','c'])

raw = np.zeros(3);

print(pLst)

pLstSorted = SortPairLst(pLst)
print(pLstSorted)

(lst1, lst2) = res = list(zip(*pLstSorted)) 
print(lst1)
print(lst2)
#print(lst3)

print(cmath.log(1j),cmath.log(-1j))

print(np.dot([1,2,3], [1,2,1]), abs(1j+1))

testLst = np.zeros( shape = (2,3,3,2) )
print( testLst[1][2] )

print( cmplxZeros( (2,3) ) )

testLst2 = np.zeros( 3 )
testLst2 += 1

print( testLst2 + 1 )

U_test = 0
U_test = ( U_test
	+ 3 
	* 22 )
	
print(U_test)

print( testFun( 1,2 ), testFunA( 3 ) )

pauliX = np.array( [[0,1],[1,0]], dtype=np.complex )
print(pauliX)

for x in range(3, 8, 2):
    print(x)
	
atu = (1,2,3)
print(atu)
print(atu + (4,5,6))
print( pauliX.shape )
print( pauliX.shape + (1,) )

arr = np.array( [[1,2,3,4], [5,6,7,8]] )
print( np.roll( arr, -1, 1 ) )