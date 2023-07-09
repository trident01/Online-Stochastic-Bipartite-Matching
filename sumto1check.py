import math
import numpy as np
from mpmath import *

mp.dps = 25
epsilon = mpf(0.11)
delta = mpf(0.18)
theta_hat = (1-epsilon)*delta/(delta + epsilon)
xLstar = 1/(1-epsilon) - delta/(delta+epsilon)


def g(x,t):
	return (1 - (1-t)*(x/(1-t) - math.floor(x/(1-t))) ) * t ** math.floor(x/(1-t))

def b(x):
	return g((1-epsilon)*x, theta_hat)

def f(x):
	return math.exp(-1*(1+delta)*x)

def fprime(x):
	return -(1+delta)*math.exp(-1*(1+delta)*x)

def bprime(x):
	return (theta_hat) ** (math.floor( (1-epsilon) * x / (1 - theta_hat))) * (epsilon - 1)

def Fprime(x):
	return b(x)*(-1 * f(1-x) + x * fprime(1-x)) + bprime(x) * (1 - x*f(1-x)) + 1 

N = 100000
print(min([ abs(Fprime(mpf(x))) for x in np.linspace(0,float(xLstar),N)[:-1]]))
print(min([ abs(Fprime(mpf(x))) for x in np.linspace(float(xLstar),1,N)[:-1]]))
