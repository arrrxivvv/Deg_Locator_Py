from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
for ind in range(20,50,5):	plt.plot( eVecLst[ind,:50,0,0] )


fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(0, ang2pi+ang_step, ang_step)
Y = X
X, Y = np.meshgrid(X, Y)
Z = quasiELst[:,:,0]
surf = ax.plot_surface(X, Y, np.real(Z), cmap=cm.coolwarm, linewidth=0, antialiased=False)
Z = quasiELst[:,:,1]
surf = ax.plot_surface(X, Y, np.real(Z), cmap=cm.coolwarm, linewidth=0, antialiased=False)
Z = quasiELst[:,:,2]
surf = ax.plot_surface(X, Y, np.real(Z), cmap=cm.coolwarm, linewidth=0, antialiased=False)
Z = quasiELst[:,:,3]
surf = ax.plot_surface(X, Y, np.real(Z), cmap=cm.coolwarm, linewidth=0, antialiased=False)

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(0, ang2pi+ang_step, ang_step)
Y = X
X, Y = np.meshgrid(X, Y)
for n in range(0,N):
   Z = np.real( quasiELst[:,:,n] )
   surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

fig = plt.figure()
plt.plot( quasiELst[20] )
fig.show()

Z = f12Lst[:,:,0]

fig = plt.figure()
ax = fig.gca(projection='3d')
Z = np.imag( f12Lst[:,:,0] )
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

# Make data.
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
