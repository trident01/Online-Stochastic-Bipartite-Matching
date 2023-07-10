import math
import numpy as np
from scipy.spatial import ConvexHull 

N = 1000 #precision parameter

def g(x,t):
	return (1 - (1-t)*(x/(1-t) - math.floor(x/(1-t))) ) * t ** math.floor(x/(1-t))

# B = beta, P = phi
def apxAtT(theta, B,P):
	if B > (1-P) ** 2:
		print("Error, beta cannot be greater than (1-phi)^2")
		return 
	if B == 0: # evaluate at the limit as B->0
		return 1 - g(P,theta)*math.exp(-1*(1-P))
	return 1-g(P, theta) *  ((1 - B/(1-P)) ** ( (1-P)*(1-P)/B))

theta = 0.5
functionPoints = []
for P in np.linspace(0,1,N):
	for B in np.linspace(0,1,N):
		if B <= (1-P) ** 2:
			functionPoints.append([B,P,apxAtT(theta, B, P)])
hull = ConvexHull(functionPoints)

def pointInHull(point):
	for eq in hull.equations:
		if np.dot(eq[:-1], point) + eq[-1] > 0:
			return False
	return True 

convexHullEvaluations = {} # for caching purposes
def evaluateAtConvexHull(theta, B, P):
	if (theta, B, P) in convexHullEvaluations:
		return convexHullEvaluations[(theta, B, P)]
	else:
		a = 0
		for i in range(15):
			if a + 2 ** (-1 * i) < apxAtT(theta, B, P) and not pointInHull([B, P, a+2 ** (-1 * i)]):
				a += 2 ** (-1 * i)
		convexHullEvaluations[(theta, B, P)] = a
		return a


def varBound(theta, S_L, S_H):
	return 1 - 0.5 * math.sqrt( S_L + S_H - S_L * (1 - theta + S_L / 2) - S_H*S_H/2)


def apxRatio(theta):
	currentMin = 1
	for S_L in np.linspace(0, theta, N):
		for S_H in np.linspace(0, 1-theta, N):
			if S_L == 0:
				print("Progress: ", 100*S_H/(1-theta))
			if varBound(theta, S_L, S_H) < currentMin:
				M = max(varBound(theta, S_L,S_H), evaluateAtConvexHull(theta, S_H*S_H/2, theta))
				if M < currentMin:
					currentMin = M
					print(currentMin)
	return currentMin

apxRatio(theta)

