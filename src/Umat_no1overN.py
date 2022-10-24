import utils

import math
import cmath
import numpy as np
from numpy import linalg as LA
import functools

def U_element_fun( N, V1, V2, tau, k1, k2, ang1, ang2 ):
	U_ans = complex(0)
	for n in range(1,N+1):
		U_ans = ( U_ans 
			+ cmath.exp(1j * 2 * math.pi * (k1-k2) * n / N ) 
			 * cmath.exp(-1j * tau * V2 * math.cos( (2*math.pi*n + ang2)/N ) ) 
			 * cmath.exp(-1j * tau * V1 * math.cos( (2*math.pi*k2 + ang1) / N )) )
	return U_ans

def U_mat_flq_raw( ang1, ang2, N, V1, V2, tau ):
	k_num = N
	U_mat = utils.cmplxZeros( (k_num, k_num) )
	for k1 in range(k_num):
		for k2 in range(k_num):
			U_mat[k1][k2] = U_element_fun( N, V1, V2, tau, k1, k2, ang1, ang2 )
	return U_mat
	
pauliX = np.array( [[0,1],[1,0]], dtype=np.complex )
pauliY = np.array( [[0,-1j],[1j,0]], dtype=np.complex )
pauliZ = np.array( [[1,0],[0,-1]], dtype=np.complex )

def U_mat_QWZ_raw( ang1, ang2, u ):
	return math.sin(ang1) * pauliX + math.sin(ang2) * pauliY + (u + math.cos(ang1) + math.cos(ang2)) * pauliZ
	
