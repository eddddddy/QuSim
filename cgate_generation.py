import numpy as np

universal_gates = [np.matrix([[0, 1], [1, 0]]), np.matrix([[1, 0], [0, -1]]), np.matrix([[0, -1j], [1j, 0]]), np.matrix([[1, 1], [1, -1]]) / np.sqrt(2), np.matrix([[1, 0], [0, 1j]]), np.matrix([[1, 0], [0, -1j]]), np.matrix([[1, 0], [0, 1/np.sqrt(2) + 1j/np.sqrt(2)]]), np.matrix([[1, 0], [0, 1/np.sqrt(2) - 1j/np.sqrt(2)]])]
rotation_gates = [np.matrix([[1, 0], [0, 1j]]), np.matrix([[1, 0], [0, -1j]]), np.matrix([[1, 0], [0, 1/np.sqrt(2) + 1j/np.sqrt(2)]]), np.matrix([[1, 0], [0, 1/np.sqrt(2) - 1j/np.sqrt(2)]])]
other_gates = [np.matrix([[0, 1], [1, 0]]), np.matrix([[1, 0], [0, -1]]), np.matrix([[0, -1j], [1j, 0]]), np.matrix([[1, 1], [1, -1]]) / np.sqrt(2)]
num_of_gates = len(universal_gates)
X = np.matrix([[0, 1], [1, 0]])
whole_representation = []
curr_A_representation = []
curr_B_representation = []
start_with_rotation = True

# determines if two matrices are "close enough" to each other (this function is required due to Python's float inprecision)
def are_close_enough(mat1, mat2):
	epsilon = 1e-6
	
	def flatten(mat):
		if len(mat) == 0:
			return []
		else:
			return mat[0] + flatten(mat[1:])
	
	return all(map(lambda x, y: np.abs(x - y) < epsilon, flatten(mat1.tolist()), flatten(mat2.tolist())))


def find_controlled_gate_AXBXC(matrix_rep):

	global whole_representation
	global curr_A_representation
	global curr_B_representation
	
	whole_representation = []
	curr_A_representation = []
	curr_B_representation = []
	
	def list_equal(l1, l2):
		if len(l1) != len(l2):
			return False
		else:
			return all(list(map(lambda x, y: x == y, l1, l2)))
			
	def iter_whole_representation():
		global whole_representation
		length = len(whole_representation)
		for i in range(length - 1, -1, -1):
			if whole_representation[i] != num_of_gates - 1:
				whole_representation[i] += 1
				for j in range(i + 1, length):
					whole_representation[j] = 0
				return
		print("Finished searching size", length)
		whole_representation = [0] * (length + 1)
	
	def modify_representations():
		global whole_representation
		global curr_A_representation
		global curr_B_representation
		if list_equal(whole_representation, curr_B_representation):
			iter_whole_representation()
			curr_A_representation = whole_representation
			curr_B_representation = []
		else:
			curr_A_representation = whole_representation[0: len(curr_A_representation) - 1]
			curr_B_representation = whole_representation[len(curr_A_representation) - 1:]
			
	def mult_all(rep):
		curr = np.matrix([[1, 0], [0, 1]])
		for i in range(len(rep)):
			curr = np.matmul(curr, universal_gates[rep[i]])
		return curr
	
	while True:
		A = mult_all(curr_A_representation)
		B = mult_all(curr_B_representation)
		C = np.linalg.inv(np.matmul(A, B))
		mul = np.matmul
		if are_close_enough(matrix_rep, mul(A, mul(X, mul(B, mul(X, C))))):
			return curr_A_representation, curr_B_representation
		modify_representations()


def find_controlled_gate_AXBX(matrix_rep):

	global whole_representation
	whole_representation = []
	
	def iter_whole_representation():
		global whole_representation
		length = len(whole_representation)
		for i in range(length - 1, -1, -1):
			if whole_representation[i] != num_of_gates - 1:
				whole_representation[i] += 1
				for j in range(i + 1, length):
					whole_representation[j] = 0
				return
		print("Finished searching size", length)
		whole_representation = [0] * (length + 1)
	
	def mult_all(rep):
		curr = np.matrix([[1, 0], [0, 1]])
		for i in range(len(rep)):
			curr = np.matmul(curr, universal_gates[rep[i]])
		return curr
	
	while True:
		A = mult_all(whole_representation)
		B = np.linalg.inv(A)
		mul = np.matmul
		if are_close_enough(matrix_rep, mul(X, mul(B, mul(X, A)))):
			return whole_representation
		iter_whole_representation()
		
		
def find_controlled_gate_AXBX_optimized(matrix_rep):
	
	global whole_representation
	global start_with_rotation
	whole_representation = []
	start_with_rotation = True
	
	def iter_whole_representation():
	
		global start_with_rotation
		
		if start_with_rotation:
			start_with_rotation = False
			return
	
		global whole_representation
		length = len(whole_representation)
		start_with_rotation = True
		
		for i in range(length - 1, -1, -1):
			if whole_representation[i] != 3:
				whole_representation[i] += 1
				for j in range(i + 1, length):
					whole_representation[j] = 0
				return
		print("Finished searching size", length)
		whole_representation = [0] * (length + 1)
		
	def mult_all_rotation(rep):
		curr = np.matrix([[1, 0], [0, 1]])
		for i in range(len(rep)):
			if i % 2 == 0:
				curr = np.matmul(curr, rotation_gates[rep[i]])
			else:
				curr = np.matmul(curr, other_gates[rep[i]])
		return curr
			
	def mult_all_other(rep):
		curr = np.matrix([[1, 0], [0, 1]])
		for i in range(len(rep)):
			if i % 2 == 0:
				curr = np.matmul(curr, other_gates[rep[i]])
			else:
				curr = np.matmul(curr, rotation_gates[rep[i]])
		return curr
			
	def mult_all(rep):
		if start_with_rotation:
			return mult_all_rotation(rep)
		else:
			return mult_all_other(rep)
	
	while True:
		A = mult_all(whole_representation)
		B = np.linalg.inv(A)
		mul = np.matmul
		if are_close_enough(matrix_rep, mul(X, mul(B, mul(X, A)))):
			return whole_representation
		iter_whole_representation()
		

def find_controlled_gate_AXB(matrix_rep):

	global whole_representation
	whole_representation = []
	
	def iter_whole_representation():
		global whole_representation
		length = len(whole_representation)
		for i in range(length - 1, -1, -1):
			if whole_representation[i] != num_of_gates - 1:
				whole_representation[i] += 1
				for j in range(i + 1, length):
					whole_representation[j] = 0
				return
		print("Finished searching size", length)
		whole_representation = [0] * (length + 1)
	
	def mult_all(rep):
		curr = np.matrix([[1, 0], [0, 1]])
		for i in range(len(rep)):
			curr = np.matmul(curr, universal_gates[rep[i]])
		return curr
	
	while True:
		A = mult_all(whole_representation)
		B = np.linalg.inv(A)
		mul = np.matmul
		if are_close_enough(matrix_rep, mul(A, mul(X, B))):
			return whole_representation
		iter_whole_representation()



















	
		
