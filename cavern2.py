# -*- coding: utf-8 -*-
import numpy as np
from box import *

class Cavern:
    #------------Model parameters-------------------
    siF = 1/3.
    siG = 1/3.
    eta= 1/3.
    L1 = 300.e3 #m
    Tinit=-2.
    Sinit=34.
    T0=1.5
    S0=34.
    #--------------------------------------------------
    
    #Constructor. Initialize the geometry. It takes typical positive depth and width which catacterize each ice shelf
    def __init__(self,W,pi,pfd,pw,pgl,Ld): 
        self.W=W
        self.pi=pi
        self.pfd=pfd
        self.pw=pw
        self.pgl = pgl
        self.Ld=Ld
        self.t=0.
        self.initGeometry()
        
    #Build the 6 box defining a cavity
    def initGeometry(self):
        pi=self.pi
        pfd = self.pfd
        pid = self.pi + self.siF * (self.pfd-self.pi)
        self.H1 = pid
        self.H2 = pfd - pid
        pw = self.pw
        pgl = self.pgl
        pwd = pw + self.siG *(pgl-pw)
        Ld = self.Ld
        Lw = self.eta*Ld
        Li = Ld - Lw
        piw = pw + Lw*(pi-pw)/Ld
        piw2 = pwd + Lw *(pid-pwd)/Ld
        W = self.W
        B1 = pid * self.L1
        B2 = self.L1 * (pfd - pid) 
        Bd = 1/2.*abs(pid-pwd)*Ld + abs(max(pgl,pfd)-max(pwd,pid))*Ld - 1/2.*abs(pfd-pgl)*Ld
        Bw = 1/2.*abs(piw-pw)*Lw + abs(max(pwd,piw2)-max(pw,piw))*Lw - 1/2.*abs(pwd-piw2)*Lw
        Bi = 1/2.*abs(piw-pi)*Li + abs(max(pid,piw2)-max(pi,piw))*Li - 1/2.*abs(pid-piw2)*Li
        
        
        self.V1=B1*W
        self.V2=B2*W
        self.Vd=Bd*W
        self.Vw=Bw*W
        self.Vi=Bi*W
        self.A1 = self.L1*W
        self.A2 = self.L1*W
        self.Ad = Ld*W
        self.Aw = Lw *W
        self.Ai = Li*W
        self.Hd=Bd/Ld
        self.Hw=Bw/Lw
        self.Hi=Bi/Li
        self.Lw = Lw
        self.Li = Li
        self.pi = pi
        self.pfd = pfd
        self.pid = pid
        self.pgl = pgl
        self.pwd = pwd
        self.pw = pw
        self.piw = piw
        self.piw2 = piw2
        
        self.box1=Box(self.Tinit,self.Sinit,self.V1)
        self.box2=Box(self.Tinit,self.Sinit,self.V2)
        self.boxd=Box(self.Tinit,self.Sinit,self.Vd)
        self.boxw=Box(self.Tinit,self.Sinit,self.Vw)
        self.boxi=Box(self.Tinit,self.Sinit,self.Vi)
      
    #a time step of size dt      
    def step(self,dt,C,Kv,F1,tau,deltaRho,T0,S0):
    
    #Constructor of unknows vector
    def buildXVector(self):
      
    #Actualize box properties with system solution
    def fillXVector(self,X):
            
    def buildAMatrix(self,C,Kv,F1,tau,deltaRho):
             
    def buildFVector(self,C,F1,tau,T0,S0):
   
    #Turbulent model
    def turbuBoundary(self,T,S,H):
       
