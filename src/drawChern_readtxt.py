import numpy as np
import matplotlib.pyplot as plt
import utils

N = 9
chernFile = "clst_n9_ang20_shftSum.csv"
vFile = "V_n9_ang20_shftSum.csv"

chernLst = np.loadtxt( chernFile, delimiter=',' )
vLst = np.loadtxt( vFile, delimiter=',' )

fig, axs = plt.subplots(N+1,sharex=True,figsize=(5,10))
chernSum = np.zeros( chernLst[:,0].size )
for n in range(0,N,1):
	axs[n].plot( vLst, chernLst[:,n] )
	axs[n].set( ylabel= ("C"+str(n)) )
	chernSum += chernLst[:,n]
# plt.plot( vLst, chernLst, alpha = 0.7 )
axs[N].plot( vLst, chernSum )
axs[N].set( ylabel = "sum" )
plt.xlabel("V")
fig.show()

fig.savefig("clst_n9_sang20_shftSum.jpg")
fig.savefig("clst_n9_ang20_shftSum.eps")

input("Press Enter to continue...")