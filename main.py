import numpy as np
from box import *
from cavern import *
import matplotlib.pyplot as plt

#Run the cavity
def runCavern(PIG,C,Kv,F1,tau,deltaRho,T0,S0,draw):
    #PIG: Cavern object
    #draw: conditional of plot
    
    volTotal=PIG.Vi+PIG.Vw+PIG.Vd
    area= 1/2.*abs(pw-pi)*Ld + abs(max(pfd,pgl)-max(pi,pw))*Ld - 1/2.*abs(pfd-pgl)*Ld
    
    print volTotal, area*W  #verify areas

    
    x=[]
    y2=[]
    y1=[]
    yd=[]
    yw=[]
    yi=[]
    y2s=[]
    y1s=[]
    yds=[]
    yws=[]
    yis=[]
    yimelt=[]
    yIWmelt=[]
    ywmelt=[]
    
    dtime = 3600*8.
    counter=0
    
    dayinsec=3600*24.
    iterday=dayinsec/dtime
    
    while(counter<30000):
        counter+=1
        A=PIG.step(dtime,C,Kv,F1,tau,deltaRho,T0,S0)
        if np.mod(counter,iterday)==0:
            dia=counter/iterday
            x.append(np.log10(dia/365.))
            y2.append(PIG.box2.temp)
            y1.append(PIG.box1.temp)
            yw.append(PIG.boxw.temp)
            yd.append(PIG.boxd.temp)
            yi.append(PIG.boxi.temp)
            y2s.append(PIG.box2.salt)
            y1s.append(PIG.box1.salt)
            yws.append(PIG.boxw.salt)
            yds.append(PIG.boxd.salt)
            yis.append(PIG.boxi.salt)
            
            ywmelt.append(PIG.boxw.melt)
            yIWmelt.append((PIG.boxw.melt*PIG.Aw+PIG.boxi.melt*PIG.Ai)/(PIG.Aw+PIG.Ai))
            yimelt.append(PIG.boxi.melt)
    
    if draw:
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        ax.plot(x,y1,'blue')
        ax.plot(x,y2,'red')
        ax.plot(x,yw,'green')
        ax.plot(x,yd,'black')
        ax.plot(x,yi,'purple')
        
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.plot(x,y1s,'blue')
        ax2.plot(x,y2s,'red')
        ax2.plot(x,yws,'green')
        ax2.plot(x,yds,'black')
        ax2.plot(x,yis,'purple')
        
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.plot(x,ywmelt,'green')
        ax3.plot(x,yIWmelt,'black')
        ax3.plot(x,yimelt,'purple')
        
        plt.show()

#Optimize the Melt Rate: Given a known melt rate, it modifies C to adjust the meltRate in the cavity
def optimMeltRate(Cav,C,Kv,F1,tau,deltaRho,T0,S0,mR):
    runCavern(Cav,C,Kv,F1,tau,deltaRho,T0,S0,False)
    mRopt = (Cav.boxw.melt*Cav.Aw+Cav.boxi.melt*Cav.Ai)/(Cav.Aw+Cav.Ai)
    print mRopt
    while abs(mR-mRopt)>(mR*0.1):
        if mR<mRopt:
            C=C*0.5*(abs(mR-mRopt)>mR)+C*0.1*(abs(mR-mRopt)<=mR)
        elif mR>mRopt:
            C=C*2*(abs(mR-mRopt)>mR)+C*1.1*(abs(mR-mRopt)<=mR)
        
        
        Cav=Cavern(Cav.W,Cav.pi,Cav.pfd,Cav.pw,Cav.pgl,Cav.Ld)
        runCavern(Cav,C,Kv,F1,tau,deltaRho,T0,S0,False)
        mRopt = (Cav.boxw.melt*Cav.Aw+Cav.boxi.melt*Cav.Ai)/(Cav.Aw+Cav.Ai)
        print mRopt
    
    print 'C=',C
    return C
    
#Given the salinity at box1, it modify F1 to obtain it
def optimSeaIceF(Cav,C,Kv,F1,tau,deltaRho,T0,S0,S1):
    runCavern(Cav,C,Kv,F1,tau,deltaRho,T0,S0,False)
    S1opt = Cav.box1.salt
    print S1opt
    while abs(S1-S1opt)>(0.006):
        if S1<S1opt:
            F1=F1*0.2*(abs(S1-S1opt)>S1*0.01)+F1*0.1*(abs(S1-S1opt)<=S1*0.01)
        elif S1>S1opt:
            F1=F1*1.5*(abs(S1-S1opt)>S1*0.01)+F1*1.2*(abs(S1-S1opt)<=S1*0.01)
        
        print F1
        Cav=Cavern(Cav.W,Cav.pi,Cav.pfd,Cav.pw,Cav.pgl,Cav.Ld)
        runCavern(Cav,C,Kv,F1,tau,deltaRho,T0,S0,False)
        S1opt = Cav.box1.salt
        print S1opt
    
    print 'F=',F1
    return F1
    
def update_line(fig,h1,data):
    h1.set_xdata(np.append(h1.get_xdata(), data[0]))
    h1.set_ydata(np.append(h1.get_ydata(), data[1]))
    fig.canvas.draw()



###   PIG   ###
W=40e3
Ld=70e3
pi=400.
pfd=1000.
pw=900.
pgl=1000.

PIG=Cavern(W,pi,pfd,pw,pgl,Ld)

C = 2e6
Kv = 1e-5
F1 = 0e-10
tau =300*3600*24.
deltaRho = 0.
T0=1.1
S0=34.68
mR=26
S1=34.05

#optimMeltRate(PIG,C,Kv,F1,tau,deltaRho,T0,S0,mR)
#optimSeaIceF(PIG,C,Kv,F1,tau,deltaRho,T0,S0,S1)
runCavern(PIG,C,Kv,F1,tau,deltaRho,T0,S0,True)



###   ROS   ####
W=600e3
Ld=800e3
pi=240.
pfd=400.
pw=300.
pgl=400.

C=6859484.192
Kv = 1e-5
F1 = 0.0000000007
tau =30*3600*24.
deltaRho = 0.
T0=-1.7
S0=34.83
mR=0.13
S1=34.72
ROS=Cavern(W,pi,pfd,pw,pgl,Ld)

#optimMeltRate(ROS,C,Kv,F1,tau,deltaRho,T0,S0,mR)
#optimSeaIceF(ROS,C,Kv,F1,tau,deltaRho,T0,S0,S1)

#runCavern(ROS,C,Kv,F1,tau,deltaRho,T0,S0,True)


