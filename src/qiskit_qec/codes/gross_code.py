
from qiskit_qec.codes.qec_code import QECCode
import numpy as np

class GrossCode(QECCode):

    def __init__(self, n=144, k=12, distance=12):
        """
        Args:
            n: The code length
            k: The number of logical qubits in the code.
            d: The distance of the code.
        """
        super().__init__()

        self.n = n
        self.k = k
        self.d = distance
        


    def cyclic_shift(self):
        S = np.roll(np.eye(self.n,dtype=int), 1, axis=0)
        return S
    
    def generate_A_B(self):
        x = np.kron(self.cyclic_shift(), np.eye(self.d, dtype=int))
        y = np.kron(np.eye(self.d, dtype=int), np.eye(self.d, dtype=int))

        #A = A1 + A2 + A3, B = B1 + B2 + B3
        #A = x^3 + y + y^2, B = x + x^2 + y^3
        


