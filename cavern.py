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
    rhoE = 920. #kg/m3
    rhoStar = 1033. #Kg/m3
    alpha = 7.5e-5 #1/K
    beta = 7.7e-4
    Tfrez = -1.88 #°C
    a = -0.057 #°C/psu
    b = 0.0832 #°C
    c =7.64e-4 #°C/dbar
    Cp = 3974. #J/Kg/K
    L = 3.34e5 #J/Kg
    gammaT = 3.e-5 #m/s
    gammaS = 2.e-6 #m/s
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
        #dt: time step
        #C, Kv,F1,tau probleme parameters (see Hellmmer paper)
        #deltaRho: difference between density to induce vertical convection (tipically 0 it's OK)
        #T0, C0: ocean conditions of temp and salinity at the entrance of the box2
        
        box1=self.box1
        box2=self.box2
        boxd=self.boxd
        boxw=self.boxw
        boxi=self.boxi
        
        rhoStar=self.rhoStar
        beta=self.beta
        alpha=self.alpha
        V1=box1.volume
        V2=box2.volume
        Vw=boxw.volume
        Vi=boxi.volume
        Vd=boxd.volume
        
        H1=self.H1
        H2=self.H2
        Hw=self.Hw
        Hi=self.Hi
        Hd=self.Hd
        
        T1=box1.temp 
        S1=box1.salt
        T2=box2.temp 
        S2=box2.salt
        Td=boxd.temp 
        Sd=boxd.salt
        Tw=boxw.temp 
        Sw=boxw.salt
        Ti=boxi.temp 
        Si=boxi.salt
    
        q = C*rhoStar*(beta*(S2-Sw) - alpha * (T2 - Tw))
        q1=q/V1
        q2=q/V2
        qw=q/Vw
        qi=q/Vi
        qd=q/Vd
        
        Sbw, Tbw, mw = self.turbuBoundary(boxw.temp, boxw.salt, self.pw)
        Sbi, Tbi, mi = self.turbuBoundary(boxi.temp, boxi.salt, (self.piw+self.pi)/2)
        
        boxi.melt = mi*3600*24*365 #m/yr
        boxw.melt = mw*3600*24*365 #m/yr
        
        
        kL = Kv/100
        kH = 1.e4*kL
        k1 = kL/H1
        k2 = kL/H2
        kd = kL/Hd
        kw = kL/Hw
        ki = kL/Hi
        
        gammaT1=self.gammaT/H1
        gammaT2=self.gammaT/H2
        gammaTd=self.gammaT/Hd
        gammaTw=self.gammaT/Hw
        gammaTi=self.gammaT/Hi
        
        gammaS1=self.gammaS/H1
        gammaS2=self.gammaS/H2
        gammaSd=self.gammaS/Hd
        gammaSw=self.gammaS/Hw
        gammaSi=self.gammaS/Hi
        
        A = self.buildAMatrix(C,Kv,F1,tau,deltaRho)
        
        F = self.buildFVector(C,F1,tau,T0,S0)
        Xinit = self.buildXVector()
        
        X=np.zeros((10,1))
        
        X = Xinit +dt*(np.dot(A,Xinit)+F)
        self.fillXVector(X)
        #print X[1], A[1,:], np.dot(A,Xinit), F
        return A
    
    #Constructor of unknows vector
    def buildXVector(self):
        
        box1=self.box1
        box2=self.box2
        boxd=self.boxd
        boxw=self.boxw
        boxi=self.boxi
        
        T1=box1.temp 
        S1=box1.salt
        T2=box2.temp 
        S2=box2.salt
        Td=boxd.temp 
        Sd=boxd.salt
        Tw=boxw.temp 
        Sw=boxw.salt
        Ti=boxi.temp 
        Si=boxi.salt
        
        
        Xinit=np.zeros((10,1))
        Xinit[0]=T1
        Xinit[5]=S1
        Xinit[1]=T2
        Xinit[6]=S2
        Xinit[2]=Tw
        Xinit[7]=Sw
        Xinit[3]=Ti
        Xinit[8]=Si
        Xinit[4]=Td
        Xinit[9]=Sd
        
        return Xinit
    
    #Actualize box properties with system solution
    def fillXVector(self,X):
        
        box1=self.box1
        box2=self.box2
        boxd=self.boxd
        boxw=self.boxw
        boxi=self.boxi
        
        box1.temp=X[0] 
        box1.salt=X[5] 
        box2.temp=X[1] 
        box2.salt=X[6] 
        boxd.temp=X[4] 
        boxd.salt=X[9] 
        boxw.temp =X[2] 
        boxw.salt=X[7] 
        boxi.temp =X[3] 
        boxi.salt=X[8] 
       
       
    def buildAMatrix(self,C,Kv,F1,tau,deltaRho):
        
        box1=self.box1
        box2=self.box2
        boxd=self.boxd
        boxw=self.boxw
        boxi=self.boxi
        
        rhoStar=self.rhoStar
        beta=self.beta
        alpha=self.alpha
        V1=box1.volume
        V2=box2.volume
        Vw=boxw.volume
        Vi=boxi.volume
        Vd=boxd.volume
        
        H1=self.H1
        H2=self.H2
        Hw=self.Hw
        Hi=self.Hi
        
        rho1=box1.dens
        rho2=box2.dens
        rhod=boxd.dens
        rhow=boxw.dens
        rhoi=boxi.dens
        S2=box2.salt
        Sw=boxw.salt
        T2=box2.temp
        Tw=boxw.temp
    
        q = C * rhoStar * (beta*(S2-Sw) - alpha * (T2 - Tw))
        #print 'q' , q,C,rhoStar,S2,Sw,T2,Tw
        q1=q/V1
        q2=q/V2
        qw=q/Vw
        qi=q/Vi
        qd=q/Vd
        
        Sbw, Tbw, mw = self.turbuBoundary(boxw.temp, boxw.salt, self.pw)
        Sbi, Tbi, mi = self.turbuBoundary(boxi.temp, boxi.salt, self.piw)
        
       
        kL = Kv/100.
        kH = 1.e4*kL
        k1 = kL/H1
        k2 = kL/H2
        kw = kL/Hw
        ki = kL/Hi
        
        gammaTw=self.gammaT/Hw
        gammaTi=self.gammaT/Hi
        
        gammaSw=self.gammaS/Hw
        gammaSi=self.gammaS/Hi
        
        mw=mw/Hw
        mi=mi/Hi
        
        if (rho1-rho2)>deltaRho:
            print 'HHHHELLLO', (rho1-rho2)
            
        if (rhow-rhod)>deltaRho:
            print 'HHHHELLLO', (rhow-rhod)
            
        if (rhoi-rhod)>deltaRho:
            print 'HHHHELLLO', (rhoi-rhod)
        
        A=np.zeros((10,10))
        
        A[(0,5),(0,5)]+= -q1 -k1*((rho1-rho2)<=deltaRho) -kH*((rho1-rho2)>deltaRho)
        A[(1,6),(1,6)]+= -q2 -k2*((rho1-rho2)<=deltaRho) -kH*((rho1-rho2)>deltaRho)
        A[(2,7),(2,7)]+= -qw -kw*((rhow-rhod)<=deltaRho) -kH*((rhow-rhod)>deltaRho)
        A[(3,8),(3,8)]+= -qi -ki*((rhoi-rhod)<=deltaRho) -kH*((rhoi-rhod)>deltaRho)
        A[(4,9),(4,9)]+= -qd -(Vi/Vd)*ki*((rhoi-rhod)<=deltaRho) -(Vw/Vd)*kw*((rhow-rhod)<=deltaRho) -(Vi/Vd)*kH*((rhoi-rhod)>deltaRho) -(Vw/Vd)*kH*((rhow-rhod)>deltaRho)
        
        A[0,0] += -1/tau
        
        A[(0,5),(3,8)] += q1
        A[(2,7),(4,9)] += qw
        A[(3,8),(2,7)] += qi
        A[(4,9),(1,6)] += qd
        
        A[(0,5),(1,6)] += k1*((rho1-rho2)<=deltaRho) +kH*((rho1-rho2)>deltaRho)
        A[(1,6),(0,5)] += k2*((rho1-rho2)<=deltaRho) +kH*((rho1-rho2)>deltaRho)
        A[(2,7),(4,9)] += kw*((rhow-rhod)<=deltaRho) +kH*((rhow-rhod)>deltaRho)
        A[(3,8),(4,9)] += ki*((rhoi-rhod)<=deltaRho) +kH*((rhoi-rhod)>deltaRho)
        A[(4,9),(3,8)] += (Vi/Vd)*(ki*((rhoi-rhod)<=deltaRho)+kH*((rhoi-rhod)>deltaRho))
        A[(4,9),(2,7)] += (Vw/Vd)*(kw*((rhow-rhod)<=deltaRho)+kH*((rhow-rhod)>deltaRho))
        
        A[2,2]  += -(gammaTw+mw)
        A[3,3]  += -(gammaTi+mi)
        A[7,7]  += -(gammaSw+mw)
        A[8,8]  += -(gammaSi+mi)
        
        
        #print -q2 -k2*((rho1-rho2)<=deltaRho) -kH*((rho1-rho2)>deltaRho), k2*((rho1-rho2)<=deltaRho) +kH*((rho1-rho2)>deltaRho),q2,k2
        
        #print 'A8', A[8,8], 'A2' , A[2,2] , 'A3' ,A[3,3], 'A7',A[7,7] 
        
        return A
        
    def buildFVector(self,C,F1,tau,T0,S0):
        
        box2=self.box2
        boxw=self.boxw
        boxi=self.boxi
        
        rhoStar=self.rhoStar
        beta=self.beta
        alpha=self.alpha
        V2=box2.volume
        Hw=self.Hw
        Hi=self.Hi
        
        S2=box2.salt
        Sw=boxw.salt
        T2=box2.temp
        Tw=boxw.temp
    
        q = C*rhoStar*(beta*(S2-Sw) - alpha * (T2 - Tw))
        q2=q/V2
        Sbw, Tbw, mw = self.turbuBoundary(boxw.temp, boxw.salt, self.pw)
        Sbi, Tbi, mi = self.turbuBoundary(boxi.temp, boxi.salt, self.piw)
        
        gammaTw=self.gammaT/Hw
        gammaTi=self.gammaT/Hi
        gammaSw=self.gammaS/Hw
        gammaSi=self.gammaS/Hi
        
        mw=mw/Hw
        mi=mi/Hi
        
        F=np.zeros((10,1))
        F[0]=self.Tfrez/tau
        F[5]=F1
        F[1]=q2*T0
        F[6]=q2*S0
        F[2]=(gammaTw+mw)*Tbw
        F[7]=(gammaSw+mw)*Sbw
        F[3]=(gammaTi+mi)*Tbi
        F[8]=(gammaSi+mi)*Sbi
        
        #print 'F8', F[8], 'F2' , F[2] , 'F3' ,F[3], 'F7',F[7] ,Sbi,Tbi
        #print 'F1', F[1]
        
        return F
        
    #Turbulent model
    def turbuBoundary(self,T,S,H):
        rhoStar=self.rhoStar
        a=self.a
        b=self.b
        c=self.c
        
        press = rhoStar * 9.8 * H #Pa
        
        press = press/1.e4 #dBar
        lamb = self.L/self.Cp
        ni = self.rhoE/rhoStar
        d = lamb * self.gammaS/self.gammaT
       
        Sb = ((T-b+c*press+d) - np.sqrt((T-b+c*press+d)**2-4*a*d*S))/(2*a)
        Tb = a*Sb + b - c*press
        mb = self.gammaT*(Tb-T)/(-ni*lamb)
        return Sb,Tb,mb
