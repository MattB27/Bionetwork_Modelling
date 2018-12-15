''' SBML file converted to Gekko Format
This model assumes a time basis of seconds, events occuring in sequential order and ignores units.
'''
import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt

m = GEKKO()
m.time = np.linspace(0,100,101)
# Variables
N = m.Var(value = 1.0)
S = m.Var(value = 0.0588235)
I1 = m.Var(value = 0.00176967)
I2 = m.Var(value = 1e-06)
R1 = m.Var(value = 0.439407)
R2 = m.Var(value = 0.0)
V = m.Var(value = 0.9)
Iv2 = m.Var(value = 0.5)
J2 = m.Var(value = 0.0)
J1 = m.Var(value = 0.0)
R = m.Var(value = 0.0)

# Parameters
l_e = m.Param(value = 50.0)
R0 = m.Param(value = 17.0)
p = m.Param(value = 1.0)
tau = m.Param(value = 0.7)
theta = m.Param(value = 0.5)
nu = m.Param(value = 0.5)
eta = m.Param(value = 0.5)
tInf = m.Param(value = 21.0)
tImm = m.Param(value = 20.0)
tImm_V = m.Param(value = 50.0)

# Intermediates
mu = m.Intermediate(1 / l_e)
beta = m.Intermediate(R0 * (gamma + mu))
gamma = m.Intermediate(365 / tInf)
sigma = m.Intermediate(1 / tImm)
sigmaV = m.Intermediate(1 / tImm_V)
strain1_frac = m.Intermediate((I1 + J1) / N)
strain2_frac = m.Intermediate((I2 + J2 + Iv2) / N)
S_frac = m.Intermediate(S / N)
V_frac = m.Intermediate(V / N)
R_1_frac = m.Intermediate((R1 + R) / N)
R_2_frac = m.Intermediate((R2 + R) / N)
R_frac = m.Intermediate(R / N)
r1 = m.Intermediate(mu * (1 - p) * N)
r2 = m.Intermediate(mu * p * N)
r3 = m.Intermediate(mu * S)
r4 = m.Intermediate(mu * V)
r5 = m.Intermediate(mu * I1)
r6 = m.Intermediate(mu * I2)
r7 = m.Intermediate(mu * Iv2)
r8 = m.Intermediate(mu * R1)
r9 = m.Intermediate(mu * R2)
r10 = m.Intermediate(mu * J1)
r11 = m.Intermediate(mu * J2)
r12 = m.Intermediate(mu * R)
r13 = m.Intermediate(beta * S * ((I1 + J1) / N))
r14 = m.Intermediate(beta * S * ((I2 + J2 + Iv2) / N))
r15 = m.Intermediate(beta * (1 - tau) * V * ((I2 + J2 + Iv2) / N))
r16 = m.Intermediate(gamma * I1)
r17 = m.Intermediate(gamma * I2)
r18 = m.Intermediate(beta * (1 - theta) * R2 * (I1 + J1) / N)
r19 = m.Intermediate(beta * (1 - theta) * R1 * (I2 + J2 + Iv2) / N)
r20 = m.Intermediate(gamma / (1 - nu) * J1)
r21 = m.Intermediate(gamma / (1 - nu) * J2)
r22 = m.Intermediate(gamma / (1 - eta) * Iv2)
r23 = m.Intermediate(sigma * R1)
r24 = m.Intermediate(sigma * R2)
r25 = m.Intermediate(sigma * R)
r26 = m.Intermediate(sigmaV * V)

# Equations
m.Equation(N.dt() == (0))
m.Equation(S.dt() == ((1.0 * r1) + (-1.0 * r3) + (-1.0 * r13) + (-1.0 * r14) + (1.0 * r23) + (1.0 * r24) + (1.0 * r25) + (1.0 * r26)))
m.Equation(I1.dt() == ((-1.0 * r5) + (1.0 * r13) + (-1.0 * r16)))
m.Equation(I2.dt() == ((-1.0 * r6) + (1.0 * r14) + (-1.0 * r17)))
m.Equation(R1.dt() == ((-1.0 * r8) + (1.0 * r16) + (-1.0 * r19) + (-1.0 * r23)))
m.Equation(R2.dt() == ((-1.0 * r9) + (1.0 * r17) + (-1.0 * r18) + (-1.0 * r24)))
m.Equation(V.dt() == ((1.0 * r2) + (-1.0 * r4) + (-1.0 * r15) + (-1.0 * r26)))
m.Equation(Iv2.dt() == ((-1.0 * r7) + (1.0 * r15) + (-1.0 * r22)))
m.Equation(J2.dt() == ((-1.0 * r11) + (1.0 * r19) + (-1.0 * r21)))
m.Equation(J1.dt() == ((-1.0 * r10) + (1.0 * r18) + (-1.0 * r20)))
m.Equation(R.dt() == ((-1.0 * r12) + (1.0 * r20) + (1.0 * r21) + (1.0 * r22) + (-1.0 * r25)))

# Global Options
m.options.CV_TYPE = 0
m.options.COLDSTART = 1
m.options.IMODE = 7
m.GUI = True


m.solve()