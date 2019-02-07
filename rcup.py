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
		

# class name is AlgebraicNumber, but this only implements numbers in Z[i, exp(i * pi / 4)]	
class AlgebraicNumber:
	
	def __init__(self, a_0, b_0, a_1, b_1):
		self.a_0 = a_0
		self.b_0 = b_0
		self.a_1 = a_1
		self.b_1 = b_1
		
	def negate(self):
		return AlgebraicNumber(-a_0, -b_0, -a_1, -b_1)
		
	def add(self, other_algebraic_number):
		return AlgebraicNumber(a_0 + other_algebraic_number.a_0,
		                       b_0 + other_algebraic_number.b_0,
							   a_1 + other_algebraic_number.a_1,
							   b_1 + other_algebraic_number.b_1)
							   
	def norm_squared(self):
		return AlgebraicNumber(a_0 ** 2 + 2 * b_0 ** 2 + a_1 ** 2 + 2 * b_1 ** 2,
		                       2 * a_0 * b_0 + 2 * a_1 * b_1, 0, 0)


'''
@param r: float

@return: boolean, float

assimilate returns an integer if r is close to one
'''
def assimilate(r):
	epsilon = 10e-5
	lower = int(np.floor(r))
	upper = int(np.ceil(r))
	if np.abs(r - lower) < epsilon:
		return True, lower
	elif np.abs(r - upper) < epsilon:
		return True, upper
	else:
		return False, r


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
		epsilons = list(set(map(lambda tuple: tuple[0], A)))
		for epsilon in epsilon:
			Delta = []
			for tuple in A:
			
	
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
TODO: write this function
'''
def min_t_count(x_prime, m, k):
	return 9000

















