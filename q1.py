import numpy as np
from scipy.optimize import fsolve
from heat import LinReg

#CONSTANTES
Tfin = 150 + 273.15 #K
mfdot = 1000 / 3600 #Kg/s
L = 30
W = 2.5
qsolar = 750 * L * W

def iter(z):
    Tfout = z[0]
    Tcp = (Tfout + Tfin) / 2
    cp = LinReg(2.040*10**3, 2.067*10**3, 473, 483, Tcp)
    F = np.empty((1))
    F[0] = mfdot*cp*(Tfout - (150 + 273.15)) - qsolar
    return F

initGuess = np.array([500])
solution = fsolve(iter, initGuess)

print("(Q1) Tfout =", solution, "[K]")
