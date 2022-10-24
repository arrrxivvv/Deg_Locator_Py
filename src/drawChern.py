import numpy as np
import matplotlib.pyplot as plt
import utils

N = 5
chernFile = "Clst_n5_V1_8_V2_5_fine.npy"
vFile = "Vlst_n5_V1_8_V2_5_fine.npy"

chernLst = np.load( chernFile )
# vLst = np.load( vFile )

cnBnd = int(np.max(chernLst))+1

# fig, axs = plt.subplots(N+1,sharex=True,figsize=(5,10))
fileTitle = "cContourFine_n_"
chernSum = np.zeros( chernLst[:,0].size )
for n in range(0,N,1):
	fig = plt.figure()
	plt.contourf( vLst[0], vLst[1], chernLst[:,:,n].T, 20, cmap='seismic', vmin = -cnBnd, vmax = cnBnd )
	plt.colorbar( ticks=range( -cnBnd, cnBnd, 1 ) )
	fig.suptitle( str(n+1) )
	plt.xlabel("V_1")
	plt.ylabel("V_2")
	fig.show()
	fig.savefig(fileTitle + str(N) + "_" + str(n) + ".jpg")
	fig.savefig(fileTitle + str(N) + "_" + str(n) + ".eps")
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