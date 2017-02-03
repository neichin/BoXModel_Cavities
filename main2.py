import numpy as np
from box import *
from cavern2 import *
import matplotlib.pyplot as plt


###   PIG   ###
W=40e3
Ld=70e3
pi=400.
pfd=1000.
pw=900.
pgl=1000.

PIG=Cavern(W,pi,pfd,pw,pgl,Ld)

#-------Comprobation of cavern geometry--------#
#---------------------------------------------------------------#
# Total volume of the cavern (area in the figure * W)-----------#
# must match the addition of three interior box volumes---------#
#---------------------------------------------------------------# 
volTotal=PIG.Vi+PIG.Vw+PIG.Vd
area = 1/2.*abs(pw-pi)*Ld + abs(max(pfd,pgl)-max(pi,pw))*Ld - 1/2.*abs(pfd-pgl)*Ld
print volTotal, area*W  #verify areas



###   ROS   ####
W=600e3
Ld=800e3
pi=240.
pfd=400.
pw=300.
pgl=400.

ROS=Cavern(W,pi,pfd,pw,pgl,Ld)

#---------------------------------------------------------------#
# Total volume of the cavern (area in the figure * W)-----------#
# must match the addition of three interior box volumes---------#
#---------------------------------------------------------------# 
volTotal=ROS.Vi+ROS.Vw+ROS.Vd
area = 1/2.*abs(pw-pi)*Ld + abs(max(pfd,pgl)-max(pi,pw))*Ld - 1/2.*abs(pfd-pgl)*Ld
print volTotal, area*W  #verify areas

