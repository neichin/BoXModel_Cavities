# projection class
class Box:
    alpha = 7.7e-4 #1/K
    beta = 7.7e-4
    rhoStar = 1033 #Kg/m3
    TStar = 0
    SStar = 34 #psu
    H = 0.
    V = 0.
    W = 0.
    def __init__(self,T,S,V):
        self.temp=T
        self.salt=S
        self.volume=V
        self.dens = self.rhoStar*(1-self.alpha*(T-self.TStar) + self.beta*(S-self.SStar))
        self.melt = 0