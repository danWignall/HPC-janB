from graphics import *
import sys
import numpy as np

#
#create image of final file
#

#iter=int(sys.argv[1])
N=1000#int(sys.argv[2])
finwin=GraphWin("box",1000,1000)
finwin.setBackground("white")

fileName=sys.argv[1]
pot=np.loadtxt("."+fileName+"Big.dat")

##rescale potential to [0,255] for colours
mx=0.0 #max absolute value of pot
for row in pot:
    mx=max(mx,row[2],-row[2])
for row in pot:
    row[2]=255*row[2]/mx


#create list of points to colour (only colour points with non zero pot)
p=[]
for i in range(0,N,int(N/1000)):
    for j in range(0,N,int(N/1000)):
        if(abs(pot[i+N*j][2])>10):
            p.append(pot[i+N*j])

#colour points
r=[]       
for i in range(0,len(p)):
    r.append(Point(p[i][0]/(1000/N),p[i][1]/(1000/N)))
    r[i].setFill(color_rgb(max(0.0,p[i][2]),0,-min(0.0,p[i][2])))
    r[i].draw(finwin)

# saves the current TKinter object in postscript format
finwin.postscript(file="image.eps")#, colormode='color')

# Convert from eps format to gif format using PIL
from PIL import Image as NewImage
finwin = NewImage.open("image.eps")
finwin.save("."+fileName+"Big.gif", "gif")

#finwin.getMouse() # Pause to view result
