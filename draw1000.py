from graphics import *
import sys
import numpy as np

#
#create image of final file
#

iter=int(sys.argv[1])
N=1000#int(sys.argv[2])
finwin=GraphWin("box",1000,1000)
finwin.setBackground("white")


        
        
pot=np.loadtxt(".bigphi.dat")

##rescale potential
mx=0.0 #max absolute value of pot, used to rescale down to [0,255] for colours
for row in pot:
    mx=max(mx,row[2],-row[2])
for row in pot:
    row[2]=255*row[2]/mx


#create list of points to colour (only colour points with non zero pot)
p=[]
for i in range(0,N,int(N/1000)):
    for j in range(0,N,int(N/1000)):
        #print pot[i+N*j]
        if(abs(pot[i+N*j][2])>10):
            p.append(pot[i+N*j])

r=[]       

for i in range(0,len(p)):
    r.append(Point(p[i][0]/(1000/N),p[i][1]/(1000/N)))
    r[i].setFill(color_rgb(max(0.0,p[i][2]),0,-min(0.0,p[i][2])))
    r[i].draw(finwin)



"""
#create rectangles to colour
r=[]
for i in range(0,1000):
    r.append([])
    for j in range(0,1000):
        if(abs(pot[i][2])<1e-6):
            continue
        r[i].append(Point(i,j))
        r[i][j].draw(finwin)

#colour points according to potential, +ve pot red, -ve pot blue
for i in range (0,N,N/1000):
    if(abs(pot[i][2])<1e-6):
        continue
    
    r[int(pot[i][0])][int(pot[i][1])].setFill(color_rgb(max(0.0,pot[i][2]),0,-min(0.0,pot[i][2])))

"""
finwin.getMouse() # Pause to view result
