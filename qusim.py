####################################### Notes #######################################
#
# This module is an (attempted) implementation of a quantum computer simulator. It
# provides several useful classes (QuBit, QuGate, QuScore) to run the simulation. 
# Note that we only use 9 gates in this implementation, which are our universal gates.
# Also note that our multi-qubit gates can only operate on adjacent qubits in a
# specific orientation. (However, we can define new gates using a composition of the
# universal gates for non-adjacent qubits or other orientations.)
#
# 
#
#####################################################################################

import numpy as np

class QuBit:

	# initializes the qubit to a classical "0" decohered state
	def __init__(self):
		self.row_vector = np.matrix([np.complex(1,0),np.complex(0,0)])
		self.col_vector = np.transpose(self.row_vector)
	
	# modifies theta (angle between Bloch vector and +Z) and phi (angle between Bloch vector and +X)
	def recalculate_phases(self):
		Pass # do this later when i actually understand whats happening

	# modifies the qubit based on application of the gate
	def apply_gate(self, gate):
		self.col_vector = np.matmul(gate.matrix, self.col_vector)
		self.row_vector = np.tranpose(self.col_vector)
		recalculate_phases()
		

# stores a name and a numpy matrix
class QuGate:

	# initializes the gate to store its name , number of inputs, and matrix representation
	# note: if name is not one corresponding to a universal gate, then the caller must
	#   must specify score_to_cast
	def __init__(self, name, score_to_cast = None):
		self.name = name
		if name == "I": 
			self.gate = np.matrix([[1, 0], [0, 1]])
			self.num_of_qubits = 1
		elif name == "X": 
			self.gate = np.matrix([[0, 1], [1, 0]])
			self.num_of_qubits = 1
		elif name == "Z": 
			self.gate = np.matrix([[1, 0], [0, -1]])
			self.num_of_qubits = 1
		elif name == "Y": 
			self.gate = np.matrix([[0, -1j], [1j, 0]])
			self.num_of_qubits = 1
		elif name == "H": 
			self.gate = np.matrix([[1, 1], [1, -1]]) / np.sqrt(2)
			self.num_of_qubits = 1
		elif name == "S": 
			self.gate = np.matrix([[1, 0], [0, 1j]])
			self.num_of_qubits = 1
		elif name == "St": 
			self.gate = np.matrix([[1, 0], [0, -1j]])
			self.num_of_qubits = 1
		elif name == "CNOT": 
			self.gate = np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
			self.num_of_qubits = 2
		elif name == "T": 
			self.gate = np.matrix([[1, 0], [0, 1/np.sqrt(2) + 1j/np.sqrt(2)]])
			self.num_of_qubits = 1
		elif name == "Tt": 
			self.gate = np.matrix([[1, 0], [0, 1/np.sqrt(2) - 1j/np.sqrt(2)]])
			self.num_of_qubits = 1
		else:
			self.gate = score_to_cast.get_unitary_matrix()
			self.num_of_qubits = score_to_cast.get_num_of_qubits()
		
	# returns the name of this gate
	def get_name(self):
		return self.name
		
	# returns the number of input qubits
	def get_num_of_qubits(self):
		return self.num_of_qubits
		
	# returns the unitary matrix corresponding to an application of this gate
	def get_matrix_rep(self):
		return self.gate	

		
class QuScore:

	# class-wide constants for the assimilation function
	epsilon = 1e-8
	root_two_inverse = 1/np.sqrt(2)
	root_three_inverse = 1/np.sqrt(3)

	# example score with 3 qubits:
	# {1: {1: (H, 1), 3: (H, 1)}, 2: {1: (CNOT, 1), 2: (CNOT, 2)}}

	# initializes the score to contain the number of qubits and the current
	#   dictionary (which is empty)
	def __init__(self, num_of_qubits):
		self.num_of_qubits = num_of_qubits
		self.score = {}
		
	# returns the number of qubits in the score
	def get_num_of_qubits(self):
		return self.num_of_qubits
		
	# returns the current score as a dictionary
	def get_score(self):
		return self.score
		
	# adds a qubit to the score, initialized to the 0 state with no gates applied
	def add_qubit(self):
		self.num_of_qubits += 1
		return self
	
	# adds a gate to the specified qubit(s); it is the caller's responsibility to
	#   ensure that the qubit exists in the score and that the gate does not
	#   conflict with any other gate
	# note: multi-qubit gates operate on the passed qubit and the immediate next ones
	def add_gate(self, gate, time_step, qubit1):
	
		if gate.get_num_of_qubits() == 1:
		
			if time_step not in self.score.keys():
				self.score[time_step] = {}
			self.score[time_step][qubit1] = (gate, 1)
			
		elif gate.get_num_of_qubits() == 2:
		
			if time_step not in self.score.keys():
				self.score[time_step] = {}
			self.score[time_step][qubit1] = (gate, 1)
			self.score[time_step][qubit1 + 1] = (gate, 2)
			
		elif gate.get_num_of_qubits() == 3:
			
			if time_step not in self.score.keys():
				self.score[time_step] = {}
			self.score[time_step][qubit1] = (gate, 1)
			self.score[time_step][qubit1 + 1] = (gate, 2)
			self.score[time_step][qubit1 + 2] = (gate, 3)
			
		return self
			
	# appends a different score to the end of this one; the score passed is unmodified
	def append_score(self, following_score):
	
		last_time_step = 0
		try:
			last_time_step = max(self.score.keys())
		except ValueError:
			pass

		for time_step in following_score.score:
			self.score[time_step + last_time_step] = following_score.score[time_step]

		return self
		
	# returns an assimilated matrix with the same entries as the input matrix
	def assimilate(self, matrix):
		
		# TODO: add more collection points
		def assimilate_float(float):
			for i in range(-18, 19):
				if np.abs(float - i/6) < self.epsilon:
					return i/6
			for i in range(-4, 5):
				if np.abs(float - i * self.root_two_inverse) < self.epsilon:
					return i * self.root_two_inverse
			for i in range(-5, 6):
				if np.abs(float - i * self.root_three_inverse) < self.epsilon:
					return i * self.root_three_inverse
			return float
			
		return np.matrix(list(map(lambda l: list(map(assimilate_float, l)), matrix.tolist())))
			
		
	# helper function which returns the unitary matrix representation of the current score
	def get_unitary_matrix(self):
	
		# returns a unitary matrix corresponding to applications of the
		#   individual gates
		def get_unitary_matrix_helper(matrices):
			if len(matrices) == 1:
				return matrices[0]
			elif matrices[0] is None:
				return get_unitary_matrix_helper(matrices[1::])
			else:
				return self.assimilate(np.kron(matrices[0], get_unitary_matrix_helper(matrices[1::])))
	
		total_time_steps= len(self.score)
		time_steps_count = 0
		current_time = 1
		current_matrix = np.matrix(np.eye(2 ** self.num_of_qubits))
		
		while time_steps_count < total_time_steps:
		
			if current_time not in self.score:
				current_time += 1
				continue
				
			operations = self.score[current_time]
			matrices = []
			
			for qubit in range(1, self.num_of_qubits + 1):
				if qubit not in operations:
					matrices.append(I.get_matrix_rep())
				else:
					gate = operations[qubit][0]
					if gate.get_num_of_qubits() == 1:
						matrices.append(gate.get_matrix_rep())
					elif gate.get_num_of_qubits() == 2:
						if operations[qubit][1] == 1:
							matrices.append(None)
						elif operations[qubit][1] == 2:
							matrices.append(gate.get_matrix_rep())
					elif gate.get_num_of_qubits() == 3:
						if operations[qubit][1] == 1:
							matrices.append(None)
						elif operations[qubit][1] == 2:
							matrices.append(None)
						elif operations[qubit][1] == 3:
							matrices.append(gate.get_matrix_rep())
			
			current_matrix = self.assimilate(np.matmul(get_unitary_matrix_helper(matrices), current_matrix))
			time_steps_count += 1
			current_time += 1
		
		return current_matrix
			
	# returns the state after the score is complete, starting with the ground state
	def get_end_state(self):
		return np.matmul(self.get_unitary_matrix(), 
		                 np.transpose(np.matrix([[1] + [0] * (2 ** self.num_of_qubits - 1)])))
	
	# returns a list corresponding to the occurence probability of each classical state
	#   at the end of the score
	def measure(self):
		return list(map(lambda x: (x[0].__abs__()) ** 2, self.get_end_state.tolist()))
		

####################################### Common Gates #######################################

# one-bit identity gate
I = QuGate("I")

#################### Universal Gates ####################

# a bit-flip; also a pi rotation around the X-axis on the Bloch sphere
X = QuGate("X")
			   
# a phase-flip; also a pi rotation around the Z-axis on the Bloch sphere
Z = QuGate("Z")
			   
# both a bit-flip and a phase-flip; also a pi rotation around the Y-axis on the Bloch sphere
Y = QuGate("Y")
			   
# the Hadamard gate, which transforms the qubit into a state where 0 and 1 have equal
#   likelihood; also a pi rotation around the X,Z axis on the Bloch sphere
H = QuGate("H")
			   
# a right-handed half phase-flip; also a right-handed pi/2 rotation around the 
#   Z-axis on the Bloch sphere
S = QuGate("S")
			   
# a left-handed half phase-flip; also a left-handed pi/2 rotation around the Z-axis 
#   on the Bloch sphere
St = QuGate("St")
				
# the (non-commutative) controlled-not gate, which flips the second qubit iff the first
#   qubit is measured to be 1
CNOT = QuGate("CNOT")
				  
# a right-handed quarter phase-flip; also a right-handed pi/4 rotation around the
#   Z-axis on the Bloch sphere
# note: this is not a Clifford gate, but is a popular choice for making the set of
#   Clifford gates along with the CNOT gate a universal set of gates
T = QuGate("T")
			   
# a left-handed quarter phase-flip; also a left-handed pi/4 rotation around the Z-axis
#   on the Bloch sphere
# note: again, this is not Clifford but is popular for universality
Tt = QuGate("Tt")

###################### Other Gates ######################

# The ultimate goal of this module is to construct quantum circuits using only the
#   universal set of gates. Thus, we avoid implementing the following gates directly
#   using their unitary matrix representations, and instead implement them as a 
#   sequence of universal gates acting on some set of qubits

# the Kronecker product of Hadamard gates on 2 qubits
H2 = QuGate("H2", QuScore(2).add_gate(H, 1, 1).add_gate(H, 1, 2))

# the Kronecker product of Hadamard gates on 3 qubits
H3 = QuGate("H3", QuScore(3).add_gate(H, 1, 1).add_gate(H, 1, 2).add_gate(H, 1, 3))

# the Kronecker product of Hadamard gates on 4 qubits
H4 = QuGate("H4", QuScore(4).add_gate(H, 1, 1).add_gate(H, 1, 2).add_gate(H, 1, 3).add_gate(H, 1, 4))

# the controlled-Z gate, which performs a phase-flip on the second qubit iff the first
#   qubit is measured to be 1
CZ = QuGate("CZ", QuScore(2).add_gate(H, 1, 2).add_gate(CNOT, 2, 1).add_gate(H, 3, 2))

# the controlled-Y gate, which performs a bit-flip and a phase-flip on the second qubit
#   iff the first qubit is measured to be 1
CY = QuGate("CY", QuScore(2).add_gate(St, 1, 2).add_gate(CNOT, 2, 1).add_gate(S, 3, 2))

# the controlled Hadamard gate, which performs a Hadamard operation on the second qubit
#   iff the first qubit is measured to be 1
CH = QuGate("CH", QuScore(2).add_gate(H, 1, 2)
                            .add_gate(S, 2, 2)
							.add_gate(CNOT, 3, 1)
                            .add_gate(S, 4, 2)
							.add_gate(H, 5, 2)
							.add_gate(St, 6, 2)
							.add_gate(H, 7, 2)
							.add_gate(Tt, 8, 2)
							.add_gate(CNOT, 9, 1)
							.add_gate(T, 10, 2)
							.add_gate(H, 11, 2)
							.add_gate(S, 12, 2)
							.add_gate(St, 13, 1)
							.add_gate(X, 13, 2))

# the reversed controlled-not gate, which flips the first qubit iff the second qubit
#   is measured to be 1
# in terms of universal gates only:
# CNOTR = QuGate("CNOTR", QuScore(2).add_gate(H, 1, 1)
#			                        .add_gate(H, 1, 2)
# 			                        .add_gate(CNOT, 2, 1)
# 			                        .add_gate(H, 3, 1)
# 			                        .add_gate(H, 3, 2))
CNOTR = QuGate("CNOTR", QuScore(2).add_gate(H2, 1, 1).add_gate(CNOT, 2, 1).add_gate(H2, 3, 1))

# the swap gate, which swaps the states of the two input qubits
# in terms of universal gates only:
# SWAP = QuGate("SWAP", QuScore(2).add_gate(CNOT, 1, 1)
#                                 .add_gate(H, 2, 1)
#                                 .add_gate(H, 2, 2)
#                                 .add_gate(CNOT, 3, 1)
#                                 .add_gate(H, 4, 1)
#                                 .add_gate(H, 4, 2)
#                                 .add_gate(CNOT, 5, 1))
SWAP = QuGate("SWAP", QuScore(2).add_gate(CNOT, 1, 1).add_gate(CNOTR, 2, 1).add_gate(CNOT, 3, 1))

# a CNOT gate which acts on the two non-adjacent qubits in a three-qubit score (i.e. 
#   it flips the third qubit iff the first qubit is measured to be 1, regardless of
#   what the second qubit is measured to be)
# in terms of universal gates only:
# CNOT13 = QuGate("CNOT13", QuScore(3).add_gate(CNOT, 1, 2)
#                                     .add_gate(H, 2, 2)
#                                     .add_gate(H, 2, 3)
#                                     .add_gate(CNOT, 3, 2)
#                                     .add_gate(H, 4, 2)
#                                     .add_gate(H, 4, 3)
#                                     .add_gate(CNOT, 5, 2)
#                                     .add_gate(CNOT, 6, 1)
#                                     .add_gate(CNOT, 7, 2)
#                                     .add_gate(H, 8, 2)
#                                     .add_gate(H, 8, 3)
#                                     .add_gate(CNOT, 9, 2)
#                                     .add_gate(H, 10, 2)
#                                     .add_gate(H, 10, 3)
#                                     .add_gate(CNOT, 11, 2)
CNOT13 = QuGate("CNOT13", QuScore(3).add_gate(SWAP, 1, 2).add_gate(CNOT, 2, 1).add_gate(SWAP, 3, 2))

# a SWAP gate which swaps the two non-adjacent qubits in a three-qubit score and leaves
#   the other one unmodified
# in terms of universal gates only:
# SWAP13 = QuGate("SWAP13", QuScore(3).add_gate(CNOT, 1, 1)
#                                     .add_gate(H, 2, 1)
#                                     .add_gate(H, 2, 2)
#                                     .add_gate(CNOT, 3, 1)
#                                     .add_gate(H, 4, 1)
#                                     .add_gate(H, 4, 2)
#                                     .add_gate(CNOT, 5, 1)
#                                     .add_gate(CNOT, 6, 2)
#                                     .add_gate(H, 7, 2)
#                                     .add_gate(H, 7, 3)
#                                     .add_gate(CNOT, 8, 2)
#                                     .add_gate(H, 9, 2)
#                                     .add_gate(H, 9, 3)
#                                     .add_gate(CNOT, 10, 2)
#                                     .add_gate(CNOT, 11, 1)
#                                     .add_gate(H, 12, 1)
#                                     .add_gate(H, 12, 2)
#                                     .add_gate(CNOT, 13, 1)
#                                     .add_gate(H, 14, 1)
#                                     .add_gate(H, 14, 2)
#                                     .add_gate(CNOT, 15, 1)
SWAP13 = QuGate("SWAP13", QuScore(3).add_gate(SWAP, 1, 1).add_gate(SWAP, 2, 2).add_gate(SWAP, 3, 1))

# the Toffoli gates, which flips the state of the third qubit in a three-qubit score iff
#   the states of the first and second qubits are both measured to be 1
TOF = QuGate("TOF", QuScore(3).add_gate(H, 1, 3).add_gate(CNOT, 2, 2).add_gate(Tt, 3, 3)
                              .add_gate(CNOT13, 4, 1).add_gate(T, 5, 3).add_gate(CNOT, 6, 2)
							  .add_gate(Tt, 7, 3).add_gate(CNOT13, 8, 1).add_gate(T, 9, 2)
							  .add_gate(T, 9, 3).add_gate(CNOT, 10, 1).add_gate(H, 10, 3)
							  .add_gate(T, 11, 1).add_gate(Tt, 11, 2).add_gate(CNOT, 12, 1))

############################################################################################

###################################### Useful Scores #######################################

# two-qubit entanglement [0.5, 0, 0, 0.5]
entangled2_0 = QuScore(2).add_gate(H, 1, 1).add_gate(CNOT, 2, 1)

# two-qubit entanglement [0, 0.5, 0.5, 0]
entangled2_1 = QuScore(2).add_gate(H, 1, 1).add_gate(X, 1, 2).add_gate(CNOT, 2, 1)

############################################################################################































