import numpy as np

import qugate
import qureg


class QuBit:

    def __init__(self, state=None):
        if state:
            self.state = np.matrix(state).I
        else:
            self.state = np.matrix([1, 0]).I

    def __repr__(self):
        zero_prob = self.state[0][0]
        one_prob = self.state[1][0]
        return f"[{zero_prob} {one_prob}]"

    def format(self, form):
        zero_prob = self.state[0][0]
        one_prob = self.state[1][0]
        if form == 'dirac':
            return f"{zero_prob}|0> + {one_prob}|1>"
        elif form == 'row':
            return f"[{zero_prob} {one_prob}]"
        elif form == 'long':
            return f"Coefficient of 0: {zero_prob}\nCoefficient of 1: {one_prob}"

class QuBits:

    def __init__(self, *qubits):
