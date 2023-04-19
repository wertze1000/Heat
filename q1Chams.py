import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


Tfin = 150+273.15
mf = 1000/3600

#Géométrie du problème

L = 30
W = 2.5

#droite du cp
m = (1.905*10**3-1.877*10**3)/(423-413)
p = 1.905*10**3-m*423


qsolar = 750*L*W
qconv = qsolar

def temp(z):
    Tfout = z[0]
    
    Tcp = (Tfin + Tfout)/2
    cp = m*Tcp + p
    print("CP", cp, "slope", m, "p", p)
    F = np.empty((1))
    
    F[0] = mf*cp*(Tfout-Tfin)-qconv
    
    return F

zGuess = np.array([500])
z = fsolve(temp, zGuess)

print(z)