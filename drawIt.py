
from graphics import *
import sys
import numpy as np

iteration=int(sys.argv[1])
N=100

#
#iterate
#


#create image
win=GraphWin("box",1000,1000)
win.setBackground("white")

#create rectangles to colour
r=[]
for i in range(0,100):
    r.append([])
    for j in range(0,100):
        r[i].append(Rectangle(Point(i*10,j*10),Point(i*10+10,j*10+10)))
        r[i][j].setWidth(0)
        r[i][j].draw(win)

#run over iterations, colouring rectangles
for itr in range(0,iteration,5):
    if(itr==0): continue
    pot=np.loadtxt(".phi"+str(itr)+"000.dat")

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
    
    #colour rectangles
    for i in range(0,len(p)):
        r[int(p[i][0])][int(p[i][1])].setFill(color_rgb(max(0.0,p[i][2]),0,-min(0.0,p[i][2])))

#r[0][0].setFill("green")
print "drawing done"
win.getMouse() # Pause to view result
win.close() 
