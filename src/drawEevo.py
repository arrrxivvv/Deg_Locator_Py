import numpy as np
import matplotlib.pyplot as plt
import utils

EevoFile = "Eevo_n5_test.npy"

EevoLst = np.load( EevoFile )
vLst = np.load( "V_n5_V_10_test.npy" )
EevoLst2 = np.load( "Eevo_n5_test_plocal.npy" )

a = 1

EevoMax = np.max( EevoLst, -1 )
EevoMax2 = np.max( EevoLst2, -1 )

fig = plt.figure()
plt.contourf( vLst[0], vLst[1], EevoLst[:,:,0,0,0], 20, cmap='seismic' )
plt.colorbar( )
plt.xlabel("V_1")
plt.ylabel("V_2")
fig.suptitle( "thermal initial state, E_0" )
fig.show()

fig = plt.figure()
plt.contourf( vLst[0], vLst[1], EevoLst[:,:,10,10,0], 20, cmap='seismic' )
plt.colorbar( )
plt.xlabel("V_1")
plt.ylabel("V_2")
fig.suptitle( "thermal initial state, E_0, pi angle" )
fig.show()

fig = plt.figure()
plt.contourf( vLst[0], vLst[1], EevoMax[:,:,0,0], 20, cmap='seismic' )
plt.colorbar( )
plt.xlabel("V_1")
plt.ylabel("V_2")
fig.suptitle( "thermal initial state, E_fin" )
fig.show()

fig = plt.figure()
plt.contourf( vLst[0], vLst[1], EevoMax[:,:,10,10], 20, cmap='seismic' )
plt.colorbar( )
plt.xlabel("V_1")
plt.ylabel("V_2")
fig.suptitle( "thermal initial state, E_fin, pi angle" )
fig.show()

fig = plt.figure()
plt.contourf( vLst[0], vLst[1], EevoMax2[:,:,0,0], 20, cmap='seismic' )
plt.colorbar( )
plt.xlabel("V_1")
plt.ylabel("V_2")
fig.suptitle( "p0 initial state, E_fin" )
fig.show()

fig = plt.figure()
plt.contourf( vLst[0], vLst[1], EevoMax2[:,:,10,10], 20, cmap='seismic' )
plt.colorbar( )
plt.xlabel("V_1")
plt.ylabel("V_2")
fig.suptitle( "p0 initial state, E_fin, pi angle" )
fig.show()

fig = plt.figure()
plt.contourf( vLst[0], vLst[1], np.var( EevoLst2, axis=-1 )[:,:,10,10], 20, cmap='seismic' )
plt.colorbar( )
plt.xlabel("V_1")
plt.ylabel("V_2")
fig.suptitle( "p0 initial state, var(E), pi angle" )
fig.show()

fig = plt.figure()
plt.contourf( vLst[0], vLst[1], np.var( EevoLst2, axis=-1 )[:,:,0,0], 20, cmap='seismic' )
plt.colorbar( )
plt.xlabel("V_1")
plt.ylabel("V_2")
fig.suptitle( "p0 initial state, var(E)" )
fig.show()

fig = plt.figure()
plt.plot( EevoLst[0,0,0,0] )
plt.plot( EevoLst[0,0,10,10] )
plt.plot( EevoLst[0,0,20,20] )
fig.show()

fig = plt.figure()
plt.plot( EevoLst[5,5,0,0] )
plt.plot( EevoLst[5,5,10,10] )
plt.plot( EevoLst[5,5,20,20] )
fig.show()

fig = plt.figure()
plt.plot( EevoLst2[0,5,0,0] )
plt.plot( EevoLst2[0,5,10,10] )
plt.plot( EevoLst2[9,5,0,0] )
plt.plot( EevoLst2[9,5,10,10] )
fig.show()