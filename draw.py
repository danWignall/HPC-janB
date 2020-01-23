from graphics import *
import sys
import numpy as np

#
#create image of final file
#

N=100#int(sys.argv[1])
iter=int(sys.argv[1])

finwin=GraphWin("box",1000,1000)
finwin.setBackground("white")


        
        
pot=np.loadtxt("phi.dat")

##rescale potential
mx=0.0 #max absolute value of pot, used to rescale down to [0,255] for colours
for row in pot:
    mx=max(mx,row[2],-row[2])
for row in pot:
    row[2]=255*row[2]/mx


#create list of points to colour (only colour points with non zero pot)
p=[]
for i in range(0,N):
    for j in range(0,N):
        #print pot[i+N*j]
        if(abs(pot[i+N*j][2])>10):
            p.append(pot[i+N*j])

r=[]       
for i in range(0,len(p)):
    r.append(Rectangle(Point(p[i][0]*10,p[i][1]*10),Point(p[i][0]*10+10,p[i][1]*10+10)))
    r[i].setWidth(0)
    r[i].setFill(color_rgb(max(0.0,p[i][2]),0,-min(0.0,p[i][2])))
    r[i].draw(finwin)

"""
#create rectangles to colour
r=[]
for i in range(0,100):
    r.append([])
    for j in range(0,100):
        r[i].append(Rectangle(Point(i*10,j*10),Point(i*10+10,j*10+10)))
        r[i][j].setWidth(0)
        r[i][j].draw(finwin)


#colour rectangles according to potential, +ve pot red, -ve pot blue
for row in pot:    
    r[int(row[0])][int(row[1])].setFill(color_rgb(max(0.0,row[2]),0,-min(0.0,row[2])))
    

"""

finwin.getMouse() # Pause to view result


                    
