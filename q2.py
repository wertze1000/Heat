import numpy as np
from q1 import iter
from scipy.optimize import fsolve
from heat import LinReg
import math

#Q1 temps
initGuess = np.array([500])
Tfout = fsolve(iter, initGuess)
Tfin = 150 + 273.15 #K (423.15)

#Problem geometry:
W = 2.5
L = 30
d = 6*10**(-2) #cm to m
e = 2*10**(-2)
Di = d #Internal diameter
Kglass = 1.4

#External flow properties:
vInf = 3 #m/s
Tair = 25 + 273.15 #K
TmAir = ((Tfin + Tfout)/2 + Tair)/2
muAir = LinReg(184.6, 159.6, 300, 250, TmAir)*10**(-7)
kAir = LinReg(26.3, 22.3, 300, 250, TmAir)*10**(-3)

''' TABLE A4 (Air: T , rho, Cp, mu * 10**7, nu*10**6, k*10**3, alpha*10**6, Pr)
250 1.3947 1.006 159.6 11.44 22.3 15.9 0.720
300 1.1614 1.007 184.6 15.89 26.3 22.5 0.707
'''

#Internal flow properties:
mfdot = 1000 / 3600 #Kg/s

#-> Internal flow temperatures:
Tmi = Tfin
Tmo = Tfout
Tm = (Tmi + Tmo)/2 #Mean internal flow temperature (K)

mui = LinReg(0.048, 0.044, 473, 483, Tm) #mu * 10**2 /!\

print("Internal flow properties:", "Tfout =", Tfout, "Tm", Tm, "mui", mui*10**-2)

#Internal flow calculations:
ReFluid = (4 * mfdot) / (np.pi * Di * mui*10**(-2))
PrFluid = LinReg(8.4, 8.9, 473, 463 , Tm)
frictionCoeff = 1/((0.790*np.log(ReFluid) - 1.64)**2) #8.21 (p474)
NuFluid = ((frictionCoeff/8)*(ReFluid - 1000)*PrFluid)/(1 + 12.7*(np.sqrt(frictionCoeff/8))*(math.pow(PrFluid, 2/3) - 1)) #Gnielinski 8.62 (All properties at Tm)
k = LinReg(116, 117, 483, 473, Tm) #k * 10**3 /!\

hi = (NuFluid * k*10**(-3)) / Di 
print("Internal calculated numbers: ", "Reynolds", ReFluid, "Prandlt", PrFluid, "f", frictionCoeff, "Nusselt", NuFluid,"k" ,k*10**(-3), "Convection coeff hi", hi)


#External flow calculations:

#LinReg for Prandlt (Air at 25C)
PrAir = LinReg(0.720, 0.707, 250, 300, 298.15)
ReAir = (vInf*(Di + 2*e)) / (muAir)

#Nusselt with 7.52 p420 (Hilpert)
#[Re 4000–40000 -> C = 0.193 , m = 0.618]
cConstant = 0.193
mConstant = 0.618

NuAir = cConstant*(ReAir**mConstant)*math.pow(PrAir, 1/3)
hoMean = (NuAir*kAir)/(Di + 2*e)

print("External calculated numbers: ", "Reynolds", ReAir, "Prandlt", PrAir, "f", "none", "Nusselt", NuAir, "convection coeff mean ho", hoMean)
#Q3 + 2 Tf(x):

cp = LinReg(2.040,2.067,473,483, Tm)*10**3  #(LinReg at Tm for cp of internal flow)

C = (750*W)/(mfdot*cp)
B = ((1/hi) + (Di)/(2*Kglass)*np.log((Di + 2*e)/(Di)) + (Di)/(Di + 2*e)*(1/hoMean))*mfdot*cp
D = (np.pi*Di) / B
E = C + (np.pi*Di)*(1 + Tair)/B
print("Differential equation: ","f' + ", D, "f = ", E) #https://www.dcode.fr/solveur-equation-differentielle f' + D*f = E
#newTfout = 488

newTfout = 755.277 - 332.127*np.exp((-0.00725613)*30)

efficiency = (newTfout - Tfin)/(Tfout - Tfin)
print("New Tfout =", newTfout, "Efficiency =", efficiency * 100, "%")