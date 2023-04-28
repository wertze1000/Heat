import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import math
import numpy as np

#Known

T_inf = 288.15
g = 9.81
d_out = 0.1
d_in = 0.06
L = 30
k_cond = 1.4
T_f = 423.15
Ai = math.pi * d_in * L

def LinReg(Ft1, Ft2, t1, t2, tEval):
    slope = (Ft2 - Ft1) / (t2 - t1)
    p = Ft2 - slope*t2
    FtEval = slope*tEval + p
    return   FtEval

#Calcul de h_conv = hi

mu = LinReg(0.080, 0.071, 423, 433, T_f)*10**(-2)
k_conv = LinReg(123, 122,423, 433, T_f)*10**(-3)
m = 1000/3600
Re_conv = (4*m)/(math.pi*d_in*mu) #NON Laminar (>= 2300 --> fully developped)
print("ReFluid", Re_conv)
Pr_conv = LinReg(12.2, 11.2, 423, 433, T_f)
print("Prandlt fluid", Pr_conv)

#Gz_conv = d_in * Re_conv * Pr_conv/L
#Nu_conv = ((3.66/math.tanh(2.264*Gz_conv**(-1/3) + 1.7*Gz_conv**(-2/3))) + 0.0499*Gz_conv*math.tanh(Gz_conv**(-1)))/(math.tanh(2.432*Pr_conv**(1/6)*Gz_conv**(-1/6)))
#Nu_conv = 3.66 + ((0.0668*Gz_conv)/(1+0.04*Gz_conv**(2/3)))

f = (0.79*np.log(Re_conv) - 1.64)**(-2)
Nu_conv = ((f/8)*(Re_conv - 1000)*Pr_conv)/(1+12.7*(np.sqrt(f/8))*(math.pow(Pr_conv,(2/3))-1))

h_conv = Nu_conv*k_conv/d_in

#Resolution

T_s  = T_inf+1 #Guess
T_m    = 0.5*(T_s + T_inf)
k  = LinReg(18.1, 33.8, 200, 400, T_m)*10**(-3)
nu = LinReg(7.590, 26.41, 200, 400, T_m)*10**(-6)
Pr = LinReg(0.737, 0.690,200,400,T_m)
beta = 1/T_inf

Ra = ((g*beta*(T_s - T_inf)*(L**3))/(nu**2))*Pr
Nu = ((0.6 + 0.387*math.pow(Ra,(1/6))) / (math.pow(1+math.pow(0.559/Pr, 9/16), 8/27)))**2
h = Nu*k/(np.pi * L * d_out)

R_conv  = 1/(h_conv*Ai)
R_cond = math.log((d_out/d_in))/(2*math.pi*L*k_cond)
R_free = 1/h
q_eq    = (T_f - T_inf)/(R_conv + R_cond + R_free)
T_s_end = T_f - q_eq*(R_conv + R_cond)
print(T_s_end)

while abs(T_s_end - T_s) > 0.01 :
    T_s  = T_s_end
    T_m     = 0.5*(T_s + T_inf)

    k  = LinReg(18.1, 33.8, 200, 400, T_m)*10**(-3)
    nu = LinReg(7.590, 26.41, 200, 400, T_m)*10**(-6)
    Pr = LinReg(0.737, 0.690, 200, 400, T_m)
    beta = 1/T_inf
    
    Ra = ((g*beta*(T_s - T_inf)*(L**3))/(nu**2))*Pr
    Nu = ( 0.6 + 0.387*Ra**(1/6) / (1 + (0.559/Pr)**(9/16)) **(8/27) )**2
    
    h = Nu*k/L
    
    R_free = 1/h
    q_eq = (T_f - T_inf)/(R_conv + R_cond + R_free)
    T_s_end = T_f - q_eq*(R_conv + R_cond)
    print(T_s_end)
    
print("The surface temperature is", T_s_end, "heat loss:", q_eq, "[W]")