import numpy as np

import qubit
import qureg


class QuGate:

    def __init__(self, op):

        self.epsilon = 1.0e-10

        if isinstance(op, list):
            if not op:
                self.op = None
                self.size = None
            elif isinstance(op[0], QuGate):
                self.op = np.eye(op[0].size)
                for gate in op:
                    self.op *= gate.op
            elif isinstance(op[0], list):
                self.op = np.matrix(op)
            else:
                self.op = None
                self.size = None
        elif isinstance(op, np.matrix):
            self.op = op

    def __repr__(self):
        return repr(self.op)

    def __add__(self, other):
        """
        Get the tensor product of this gate and the other gate.
        This corresponds to an application of the two gates across
        adjacent qubits.
        :param other: QuGate
        :return: QuGate
        """
        result = np.kron(self.op, other.op)
        real_part, imag_part = result.real, result.imag
        real_part.flags.writeable = True
        imag_part.flags.writeable = True
        real_part[abs(real_part) < self.epsilon] = 0
        imag_part[abs(imag_part) < self.epsilon] = 0
        return QuGate(real_part + imag_part)

    def __mul__(self, other):
        """
        Get the product of this gate and the other gate.
        This corresponds to an application of the two gates across
        the same qubit.
        :param other: QuGate
        :return: QuGate
        """
        result = self.op * other.op
        real_part, imag_part = result.real, result.imag
        real_part.flags.writeable = True
        imag_part.flags.writeable = True
        real_part[abs(real_part) < self.epsilon] = 0
        imag_part[abs(imag_part) < self.epsilon] = 0
        return QuGate(real_part + imag_part)


__ops = [
    np.matrix([[1, 0], [0, 1]]),
    np.matrix([[0, 1], [1, 0]]),
    np.matrix([[0, -1j], [1j, 0]]),
    np.matrix([[1, 0], [0, -1]]),
    np.matrix([[1, 1], [1, -1]]) / np.sqrt(2),
    np.matrix([[1, 0], [0, 1j]]),
    np.matrix([[1, 0], [0, -1j]]),
    np.matrix([[1, 0], [0, np.exp(np.complex(0, np.pi / 4))]]),
    np.matrix([[1, 0], [0, np.exp(np.complex(0, -np.pi / 4))]]),
    np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
]

I, X, Y, Z, H, S, St, T, Tt, CNOT = [QuGate(__op) for __op in __ops]
