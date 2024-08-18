import sympy as smp
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import vpython

#non-dynamic variables
t, g = smp.symbols('t, g')
m1, m2 = smp.symbols('m1, m2')
L1, L2 = smp.symbols('L1, L2')

theta1, theta2 = smp.symbols('theta1, theta2', cls = smp.Function)

theta1 = theta1(t)
theta2 = theta2(t)

the1_d = smp.diff(theta1, t)
the2_d = smp.diff(theta2, t)
the1_dd = smp.diff(the1_d, t)
the2_dd = smp.diff(the2_d, t)


x1 = L1*smp.sin(theta1)
y1 = -L1*smp.cos(theta1)
x2 = L1*smp.sin(theta1)+L2*smp.sin(theta2)
y2 = -L1*smp.cos(theta1)-L2*smp.cos(theta2)


x1_f = smp.lambdify((theta1, L1), x1)
y1_f = smp.lambdify((theta1, L1), y1)
x2_f = smp.lambdify((theta1, theta2, L1, L2), x2)
y2_f = smp.lambdify((theta1, theta2, L1, L2), y2)


#Kinetic energy
K1 = 1/2 * m1 * (smp.diff(x1, t)**2 + smp.diff(y1, t)**2)
K2 = 1/2 * m2 * (smp.diff(x2, t)**2 + smp.diff(y2, t)**2)
K = K1 + K2

#Potential energy
U1 = m1*g*y1
U2 = m2*g*y2
U = U1+U2

#Lagrangian
L = K - U

#Lagrangian equations of motion
LE1 = smp.diff(L, theta1) - smp.diff(smp.diff(L, the1_d), t).simplify()
LE2 = smp.diff(L, theta2) - smp.diff(smp.diff(L, the2_d), t).simplify()

#Solving equations of motion with respect to the1_dd and the2_dd
sols = smp.solve([LE1, LE2], (the1_dd, the2_dd), simplify = False, rational = False)

#Re-parameterizing second order ODE in terms of z = dtheta/dt and turning it into numerical system
dz1dt_f = smp.lambdify((t, g, m1, m2, L1, L2, theta1, theta2, the1_d, the2_d), sols[the1_dd])
dz2dt_f = smp.lambdify((t, g, m1, m2, L1, L2, theta1, theta2, the1_d, the2_d), sols[the2_dd])
dthe1dt_f = smp.lambdify(the1_d, the1_d)
dthe2dt_f = smp.lambdify(the2_d, the2_d)

#Defining function based off the vector S and time t
#Note the use of S
def dSdt(S, t, g, m1, m2, L1, L2):
    theta1, z1, theta2, z2 = S          #Note variables are unpacked from the S list
    return [dthe1dt_f(z1), dz1dt_f(t, g, m1, m2, L1, L2, theta1, theta2, z1, z2), dthe2dt_f(z2), dz2dt_f(t, g, m1, m2, L1, L2, theta1, theta2, z1, z2)]

#obtaining time values for time series and setting constants
t = np.linspace(0, 40, 1001)
g, m1, m2, L1, L2 = 9.81, 1, 1, 2, 2


#numerically solving system and putting setting ans to list of points
#note the use of the dSdt function as a parameter, the y0 argument is passed as the first S argument into dSdt
ans = odeint(dSdt, y0 = [1, -3, -1, 5], t = t, args = (g, m1, m2, L1, L2))

#obtain the values for theta1(t) and theta2(t)
theta1 = ans.T[0]
theta2 = ans.T[2]


mass1 = vpython.sphere(color = vpython.color.red, radius = 0.3, make_trail = True)
mass2 = vpython.sphere(color = vpython.color.green, radius = 0.3, make_trail = True)

rod1 = vpython.cylinder(pos = vpython.vector(0,0,0), axis = vpython.vector(0,0,0), radius = 0.05)
rod2 = vpython.cylinder(pos = vpython.vector(0,0,0), axis = vpython.vector(0,0,0), radius = 0.05)


def getx1y1x2y2(t, theta1, theta2, L1, L2):
    return (x1_f(theta1, L1), y1_f(theta1, L1), x2_f(theta1, theta2, L1, L2), y2_f(theta1, theta2, L1, L2))

x1, y1, x2, y2, = getx1y1x2y2(t, theta1, theta2, L1, L2)

i = 0
while t[i] <= 39.9:
    vpython.rate(25)
    mass1.pos= vpython.vector(x1[i], y1[i], 0)
    mass2.pos = vpython.vector(x2[i], y2[i], 0)
    rod1.axis = vpython.vector(x1[i], y1[i], 0)
    rod2.pos = vpython.vector(x1[i], y1[i], 0)
    rod2.axis = vpython.vector(x2[i]-x1[i], y2[i]-y1[i], 0)
    i += 1







 





