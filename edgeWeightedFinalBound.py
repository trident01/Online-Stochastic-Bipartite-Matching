import math
import numpy as np
from mpmath import *

mp.dps = 25

def g(x,t):
	return (1 - (1-t)*(x/(1-t) - math.floor(x/(1-t))) ) * t ** math.floor(x/(1-t))

def bound(x,epsilon,delta,t):
	return 1- g(x*(1-epsilon),t*(1-epsilon))*math.exp(-(1+delta)*(1-x))

N = 1000000
epsilon = mpf(0.11)
delta = mpf(0.18)
t = delta/(delta + epsilon)
print(min([bound(mpf(x),epsilon, delta,t) for x in np.linspace(0,1,N)]))
