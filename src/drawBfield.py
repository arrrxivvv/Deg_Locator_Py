import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import utils

N = 2
chernFile = "Clst_n5_V1_8_V2_5_fine.npy"
vFile = "Vlst_n5_V1_8_V2_5_fine.npy"

fLst = np.load( "Bfield_sin3_ang50.npy" )
cLst = np.load( "cLst_sin3_ang50.npy" )
# fLst = np.load( "Bfield_QWZ_u1.npy" )
# cLst = np.load( "cLst_QWZ_u1.npy" )
# vLst = np.load( vFile )

# cnBnd = int(np.max(chernLst))+1

# fig, axs = plt.subplots(N+1,sharex=True,figsize=(5,10))
fileTitle = "Bfield_3sin_ang50"

fig = plt.figure()
ax = fig.gca(projection='3d')

f1 = fLst[0,:,:,:,0]
f2 = fLst[1,:,:,:,0]
f3 = fLst[2,:,:,:,0]

x = np.arange( f1.shape[0] )
y = np.arange( f1.shape[1] )
z = np.arange( f1.shape[2] )

X,Y,Z = np.meshgrid(x,y,z)

# ax.quiver(X,Y,Z, f1, f2, f3, length=0.1, normalize=True)

ang_divide = 50
ang_half = int( ang_divide/2 )

ang2pi = 2 * math.pi
ang_ln = ang_divide + 1
ang_lst = np.linspace( 0, ang2pi, ang_ln )
ang_step = ang_lst[1]

ang1grd, ang2grd, ang3grd = np.meshgrid( ang_lst, ang_lst, ang_lst )

rvec = np.stack( [ang1grd-math.pi, ang2grd-math.pi, ang3grd-math.pi] )

rsq = np.sum( np.power( rvec,2 ), axis=0 )

pole_strength = 1/2

fieldvec_pi = pole_strength * rvec / np.power( np.sum( np.power( rvec,2 ), axis=0 ), 3/2 )

fieldvec = np.zeros( fieldvec_pi.shape )

for x_sh in range(0,2):
	for y_sh in range(0,2):
		for z_sh in range(0,2):
			fieldvec_tmp = np.roll( fieldvec_pi, (x_sh*ang_half, y_sh*ang_half, z_sh*ang_half), axis=(1,2,3) )
			fieldvec += (-1)**(x_sh + y_sh + z_sh) * fieldvec_tmp

fig = plt.figure()
plt.contourf( fieldvec[2,:,:,10], 20, cmap='seismic' )
plt.colorbar( )
fig.show()

dim=0
ang3=10
fig = plt.figure()
plt.contourf( fLst[dim][:,:,ang3,0] / (ang_step**2), 20, cmap='seismic' )
plt.colorbar( )
fig.show()
# fig.suptitle( str(ang3+1) )

fig = plt.figure()
plt.plot( fLst[0][:,ang_half,ang3,0] / (ang_step**2) )
plt.plot( fieldvec[2][:,ang_half,ang3] / 4 )
fig.show()

# fig.show()
# chernSum = np.zeros( chernLst[:,0].size )
for dim in range(0,3):
	for ang3 in range(0,50,10):
		fig = plt.figure()
		# cnBnd = int( np.max(fLst[n][1][:,:,5]) )
		plt.contourf( fLst[dim][:,:,ang3,0], 20, cmap='seismic' )
		plt.colorbar( )
		fig.suptitle( str(ang3+1) )
		plt.xlabel("ang1")
		plt.ylabel("ang2")
		fig.show()
		subfilename = fileTitle + "_dim_" + str(dim) + "_ang3_" + str(ang3)
		fig.savefig(subfilename + ".jpg")
		fig.savefig(subfilename + ".eps")
		# str(n) + "_step_0.1_" + "trueVal" + 

fig, axs = plt.subplots(N+1,sharex=True,figsize=(5,10))
indCenter = int( chernLst.shape[0]/2 )+1
chernSum = np.zeros( chernLst[indCenter,:,0].size )
for n in range(0,N,1):
	axs[n].plot( vLst[1], chernLst[indCenter,:,n] )
	axs[n].set( ylabel= ("C"+str(n)) )
	chernSum += chernLst[indCenter,:,n]
# plt.plot( vLst, chernLst, alpha = 0.7 )
axs[N].plot( vLst[1], chernSum )
axs[N].set( ylabel = "sum" )
plt.xlabel("V2")
fig.show()

fig = plt.figure()
plt.contourf( vLst[0], vLst[1], np.sum( chernLst, axis=-1 ) .T, 20, cmap='seismic', vmin = -cnBnd, vmax = cnBnd )
plt.colorbar( ticks=range( -cnBnd, cnBnd, 1 ) )
fig.suptitle( str(n+1) )
plt.xlabel("V_1")
plt.ylabel("V_2")
fig.show()
# fig.savefig(fileTitle + str(N) + "_" + str(n) + ".jpg")
# fig.savefig(fileTitle + str(N) + "_" + str(n) + ".eps")

# fig.savefig("clst_n9_sang20_shftSum.jpg")
# fig.savefig("clst_n9_ang20_shftSum.eps")

input("Press Enter to continue...")