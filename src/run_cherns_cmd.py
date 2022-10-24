import numpy as np
import subprocess

N=2

python = "python"
scriptname = "Chern.py"

file_raw = "raw.csv"
file_trunc = "trunc.csv"
file_clean = "clean.csv"

V_step = 0.5
V_min = 0 + V_step
V_max = 2

print( np.arange( V_min, V_max, V_step ) )

np.savetxt( "V_list.csv", np.arange( V_min, V_max, V_step ), delimiter=',',fmt='%f' )

for V in np.arange( V_min, V_max, V_step ):
	V1 = V
	V2 = V
	command = python + " " + scriptname + " -N " + str(N) + " -V1 " + str(V1) + " -V2 " + str(V2) + " -o " + " " + file_raw + " " + file_trunc + " " + file_clean
	print( command )
	runChern = subprocess.call( command , shell=True)