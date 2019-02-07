import numpy as np

import qubit
import qugate


class QuReg:

    def __init__(self, size):
        self.size = size
        self.score = {}
        self.next_time = 1

    def __call__(self, arg):
        if isinstance(arg, np.matrix):
            return self.get_op * arg

    def add_gate(self, gate, bit=1, time=None):

        if not time:
            time = self.next_time

        num_qubits = int(np.log(gate.op.shape[0]) / np.log(2) + 1 / 2)
        time_ops = self.score[time] if time in self.score else {}

        for i in range(num_qubits):
            time_ops[bit + i] = (gate, i)

        self.score[time] = time_ops
        self.next_time = max(self.next_time, time + 1)

        print(num_qubits, ",", time_ops.keys(), ",", self.score.keys())

        return self

    def get_op(self):
        ops = [step[1] for step in sorted([[time, ops] for time, ops in self.score.items()])[::-1]]
        return np.prod([np.sum([qugate.I if bit not in op else op[bit][0] for bit in range(1, self.size + 1) if (bit in op and not op[bit][1]) or bit not in op]) for op in ops])
