import numpy as np
from scipy.optimize import fsolve
import math

#F(F -> evaluated at t1 or t2 to find F at tEval)
def LinReg(Ft1, Ft2, t1, t2, tEval):
    slope = (Ft2 - Ft1)/(t2 - t1)
    p = Ft2 - slope*t2
    Fteval = slope*tEval + p
    return Fteval       
    

#Constantes 

g = 9.81
mf = 1000/3600
pi = np.pi
kgl = 1.4


#Géométrie du problème

L = 30
W = 2.5 
d = 0.06
e = 0.02
D = d + 2*e
Ai = np.pi * d * L
Ao = D*np.pi*L
#Températures connues

Tfin = 150+273.15
Tinf_night = 15+273.15



def temp(z):
    Ts = z[0]
    Tm = (Ts+Tinf_night)/2
    
    nu_air = LinReg(7.590*10**(-6), 26.41*10**(-6),200,400,Tm)

    Pr_air = LinReg(0.737, 0.690,200,400,Tm)
    beta = 1/Tm
    RaD = ((g*beta*(Ts - Tinf_night)*(D**3))/(nu_air**2))*Pr_air
    print("RaD", RaD)
    


    NuD = (0.6 + 0.387*RaD**(1/6) / (1 + (0.559/Pr_air)**(9/16)) **(8/27) )**2
    kair = LinReg(18.1*10**(-3), 33.8*10**(-3),200,400,(Ts+Tinf_night)/2)
    ho = (kair*NuD)/D
    print("ho", ho)
    #ho = -ho
    Rconv_o = d/(ho*D)
    
    #Internal flow calcul de hi

    mui = LinReg(0.08*10**(-2), 0.071*10**(-2), 423, 433, Tfin)

    Red = 4*mf/(pi*d*mui)
    Pri = LinReg(8.4, 8.9, 473, 463, Tfin)
    f = (0.79*np.log(Red) - 1.64)**(-2)
    Nud = ((f/8)*(Red - 1000)*Pri)/(1+12.7*(f/8)**(1/2)*(Pri**(2/3)-1))
    kflu = LinReg(123*10**(-3), 122*10**(-3), 423, 433, Tfin)
    hi = (Nud*kflu)/d
    print("Reynolds",Red, hi)
    Rtoti = 1/(hi*Ai) + np.log((D/d))/(2*np.pi*L*kgl) + d/(ho*D)
    
    q = (Tfin - Tinf_night)/Rtoti
    qloss = (Tfin - Tinf_night)/Rtoti
    print("qloss",qloss)
    F = np.empty(1)
    
    F[0] = (Ts-Tinf_night)/(Rconv_o) - q
   
    
    return F

zGuess = np.array([500])
Ts = fsolve(temp, zGuess)

print(Ts)