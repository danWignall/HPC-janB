import numpy as np

N=1000000
M=950000


pot1=np.loadtxt(".big.phi1000000.dat")
pot2=np.loadtxt(".big.phi950000.dat")

big=[]
for i in range(0,len(pot2)):
    big.append(pot2[i][2]-pot1[i][2])
    print big[i]
    
print "----------"
print max(big)
