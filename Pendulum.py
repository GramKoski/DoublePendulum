import sympy as smp
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer

import Animation
import Lyapunov


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
z1_d,z2_d  = sols[the1_dd].simplify(), sols[the2_dd].simplify()

x = smp.Matrix([[theta1], [the1_d], [theta2], [the2_d]])
f = smp.Matrix([[the1_d], [z1_d], [the2_d], [z2_d]])
#print(f)
jacobian = f.jacobian(x)

print(jacobian)
    
jacobian_f = smp.lambdify((t, g, m1, m2, L1, L2, theta1, the1_d, theta2, the2_d), jacobian)

#Re-parameterizing second order ODE in terms of z = dtheta/dt and turning it into numerical system
dz1dt_f = smp.lambdify((t, g, m1, m2, L1, L2, theta1, theta2, the1_d, the2_d), sols[the1_dd])
dz2dt_f = smp.lambdify((t, g, m1, m2, L1, L2, theta1, theta2, the1_d, the2_d), sols[the2_dd])
dthe1dt_f = smp.lambdify(the1_d, the1_d)
dthe2dt_f = smp.lambdify(the2_d, the2_d)

#Defining function based off the vector S and time t
#Note the use of S
def dSdt(t, S, g, m1, m2, L1, L2):
    theta1, z1, theta2, z2 = S          #Note variables are unpacked from the S list
    return [dthe1dt_f(z1), dz1dt_f(t, g, m1, m2, L1, L2, theta1, theta2, z1, z2), dthe2dt_f(z2), dz2dt_f(t, g, m1, m2, L1, L2, theta1, theta2, z1, z2)]

#obtaining time values for time series and setting constants
t = np.linspace(0, 40, 1001)
g, m1, m2, L1, L2 = 9.81, 1, 1, 2, 2


#print(lin @ np.identity(4))
#q, r = np.linalg.qr(lin)
S0 = [np.pi/2, 0, np.pi, 0]

#numerically solving system and putting setting ans to list of points
#note the use of the dSdt function as a parameter, the y0 argument is passed as the first S argument into dSdt
ans = solve_ivp(dSdt, t_span = [0, 40], y0 = S0, t_eval = t, args = (g, m1, m2, L1, L2)).y
theta1 = ans[0]
theta2 = ans[2]

lyap = Lyapunov.getLyap(dSdt, jacobian_f, t, S0, (g, m1, m2, L1, L2))

def getx1y1x2y2(t, theta1, theta2, L1, L2):
    return (x1_f(theta1, L1), y1_f(theta1, L1), x2_f(theta1, theta2, L1, L2), y2_f(theta1, theta2, L1, L2))

x1, y1, x2, y2, = getx1y1x2y2(t, theta1, theta2, L1, L2)

Animation.animate(t, x1, y1, x2, y2)