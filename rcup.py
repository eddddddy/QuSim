'''
This module implements the RCUP solution-searching algorithm described in https://arxiv.org/pdf/1212.6964.pdf
'''

import numpy as np

class ClosedInterval:
	
	def __init__(self, lower, upper):
		self.lower = lower
		self.upper = upper
	
	def contains(self, value):
		return lower <= value and value <= upper
		
	def intersect(self, other):
		return ClosedInterval(max(lower, other.lower), min(upper, other.upper))
		
	def intersect_is_empty(self, other):
		return other.upper < lower and other.lower > upper


'''
@param r: float

@return: float

assimilate returns an integer if r is close to one
'''
def assimilate(r):
	epsilon = 10e-5
	lower = int(np.floor(r))
	upper = int(np.ceil(r))
	if np.abs(r - lower) < epsilon:
		return lower
	elif np.abs(r - upper) < epsilon:
		return upper
	else:
		return r


'''
@param n: int
@param theta: float
@param delta: float

@return: float, list (or array-like object)

RCUP_algorithm solves the Restricted Closest Unitary Problem with
T-count n, angle theta, and threshold delta

Requires: n >= 4, delta <= 1/2
'''
def RCUP_algorithm(n, theta, delta):
	m = int(np.floor((n + 1) / 2)) + 2
	L_re_0 = find_halves(np.cos(theta), m, delta)
	L_re_1 = find_halves(np.cos(theta - np.pi / 8), m, delta)
	L_im_0 = find_halves(np.sin(theta), m, delta)
	L_im_1 = find_halves(np.sin(theta - np.pi / 8), m, delta)
	I = ClosedInterval(0, alpha) # TODO: find what alpha is
	while not I.intersect_is_empty(ClosedInterval(0, delta)):
		A = []
		for tuple1 in L_re_0:
			for tuple2 in L_im_0:
				if I.intersect(ClosedInterval(0, delta)).contains(tuple1[0] + tuple2[0]):
					A.append((tuple1[0] + tuple2[0], tuple1[1], tuple1[2], tuple2[1], tuple2[2], 0))
		for tuple1 in L_re_1:
			for tuple2 in L_im_1:
				if I.intersect(ClosedInterval(0, delta)).contains(tuple1[0] + tuple2[0]):
					A.append((tuple1[0] + tuple2[0], tuple1[1], tuple1[2], tuple2[1], tuple2[2], 1))
		A = sorted(A)
	
'''
@param alpha: float
@param m: int
@param delta: float

@return: list (or array-like object)

find_halves finds all numbers of the form a + sqrt(b) that satisfy
|alpha * sqrt(2 ^ m) - (a + b * sqrt(2))| <= delta * sqrt(2 ^ m)
and a ^ 2 + 2 * b ^ 2 <= 2 ^ m
'''
def find_halves(alpha, m, delta):
	W = alpha * np.sqrt(np.power(2, -m))
	epsilon = delta * np.sqrt(np.power(2, m))
	b = int(np.floor(-np.sqrt(np.power(2, m))))
	v = alpha * np.sqrt(np.power(2, m)) - b * np.sqrt(2)
	R = []
	while b <= int(np.ceil(np.power(2, m))):
		a_min = int(np.ceil(v - epsilon))
		a_max = int(np.floor(v + epsilon))
		for a in range(a_min, a_max + 1):
			if a ** 2 + 2 * b ** 2 <= np.power(2, m):
				R.append(((v - a) * W, a, b))
		b += 1;
		v -= np.sqrt(2)
	return sorted(R)
	
'''
'''
def min_t_count(x_prime, m, k):
	c = np.power(2, m) - np.abs(x_prime) ** 2

















