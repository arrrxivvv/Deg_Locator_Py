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
			+ 1/N * cmath.exp(1j * 2 * math.pi * (k1-k2) * n / N ) 
			 * cmath.exp(-1j * tau * V2 * math.cos( (2*math.pi*n + ang2)/N ) ) 
			 * cmath.exp(-1j * tau * V1 * math.cos( (2*math.pi*k2 + ang1) / N )) )
	return U_ans
	
def U_flq_elem( N, tau, V1, V2, ang1, ang2, k1, k2 ):
	U_ans = utils.cmplxZeros( V1.shape )
	for n in range(1,N+1):
		U_ans = ( U_ans 
			+ 1/N * np.exp(1j * 2 * math.pi * (k1-k2) * n / N ) 
			 * np.exp(-1j * tau * V2 * np.cos( (2*math.pi*n + ang2)/N ) ) 
			 * np.exp(-1j * tau * V1 * np.cos( (2*math.pi*k2 + ang1) / N )) )
	return U_ans	
	
def U_mat_vec( U_elem_fun, N, tau, V1lst, V2lst, ang_divide ):
	k1lst = np.arange(N)
	k2lst = np.arange(N)
	ang2pi = 2 * math.pi
	ang_ln = ang_divide + 1
	ang_lst = np.linspace( 0, ang2pi, ang_ln )
	V1, V2, ang1, ang2, k1, k2 = np.meshgrid( V1lst, V2lst, ang_lst, ang_lst, k1lst, k2lst )
	U_ans = U_elem_fun( N, tau, V1, V2, ang1, ang2, k1, k2 )
	return U_ans

def U_mat_flq_raw( ang1, ang2, N, V1, V2, tau ):
	k_num = N
	U_mat = utils.cmplxZeros( (k_num, k_num) )
	for k1 in range(k_num):
		for k2 in range(k_num):
			U_mat[k1][k2] = U_element_fun( N, V1, V2, tau, k1, k2, ang1, ang2 )
	return U_mat
	
def U_mat_flq_raw_vec( ang1, ang2, N, V1, V2, tau ):
	k_num = N
	klst = np.arange(knum)
	k1,k2 = np.meshgrid(klst,klst)
	U_mat = U_element_fun( N,V1,V2,tau,k1,k2,ang1,ang2 )
		
pauliX = np.array( [[0,1],[1,0]], dtype=np.complex )
pauliY = np.array( [[0,-1j],[1j,0]], dtype=np.complex )
pauliZ = np.array( [[1,0],[0,-1]], dtype=np.complex )
pauliMat = np.stack( [pauliX,pauliY,pauliZ] )

def U_mat_QWZ_raw( ang1, ang2, u ):
	return math.sin(ang1) * pauliX + math.sin(ang2) * pauliY + (u + math.cos(ang1) + math.cos(ang2)) * pauliZ

def U_mat_QWZ_vec( ang1, ang2, u ):
	ang1Lst, ang2Lst, uLst = np.meshgrid( ang1, ang2, u )
	numParam = 3
	param = [0]*3
	param[0] = np.sin(ang1Lst)
	param[1] = np.sin(ang2Lst)
	param[2] = uLst + np.cos(ang1Lst) + np.cos(ang2Lst)
	Umat = 0
	for it in range(0, numParam):
		Umat += np.multiply.outer( param[it], pauliMat[it] )
	return Umat
	
def U_mat_3sin( ang1, ang2, ang3 ):
	angLst = np.meshgrid(ang1,ang2,ang3)
	Umat = 0
	numAng = 3
	for it in range(0, numAng):
		Umat += np.multiply.outer( np.sin(angLst[it]), pauliMat[it] )
	return Umat
	
