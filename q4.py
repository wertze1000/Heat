import numpy as np
from heat import LinReg
import math

TinfNight = 15 + 273.15 #K
Tfin = 150 + 273.15     #K (423.15)
mfdot = 1000 / 3600     #Kg/s
#Hypothesis Tfin ~ Tfout
kglass = 1.4            #W/m*K
Di = 6*10**(-2)         #m
e = 2*10**(-2)

#hi calculation:
Tm = Tfin #(HYPOTHESIS !)

mui = LinReg(0.071, 0.080, 433, 423, Tm)*10**(-2)
PrandltInternal = LinReg(11.2, 12.2, 433, 423, Tm)
kInternal = LinReg(122, 123, 433, 423, Tm)*10**(-3)
ReynoldsInternal = (4*mfdot)/(np.pi*Di*mui)

frictionCoeff = 1/((0.790*np.log(ReynoldsInternal) - 1.64)**2) #8.21 (p474)
NusseltInternal = ((frictionCoeff/8)*(ReynoldsInternal - 1000)*PrandltInternal)/(1 + 12.7*(np.sqrt(frictionCoeff/8))*(math.pow(PrandltInternal, 2/3) - 1))

hi = (NusseltInternal*kInternal) / Di
print("Internal flow:", " Prandlt: ", PrandltInternal, " Reynolds: ", ReynoldsInternal, " fCoeff: ", frictionCoeff, " Nusselt: ", NusseltInternal, "hi", hi)
ho = "placeholder"
#ui = 1/(1/hi + (Di)/(2*kglass)*np.log((Di + 2*e)/(Di)) + (Di)/(Di + 2*e)*(1/ho))