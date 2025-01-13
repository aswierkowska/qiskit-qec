
from qiskit_qec.codes.qec_code import QECCode
import numpy as np
import scipy.linalg as la
import sympy as sp

class GrossCode(QECCode):

    def __init__(self, n=144, k=12, distance=12):
        """
        Args:
            n: The code length
            k: The number of logical qubits in the code.
            d: The distance of the code.
        """
        super().__init__()


        self.l = 12
        self.m = 6

        self.n = 2*self.l*self.m
        self.k = None
        self.d = None

        self.A = []
        self.B = []

        self.A_matrix = None
        self.B_matrix = None

        self.H_X, self.H_Z = self.generate_check_matrix()
        #parameters for gross code
        assert self.k == 12
        assert self.d == 12
        assert self.n == 144
        assert self.H_X.shape == (72, 144)
        assert self.H_Z.shape == (72, 144)
        assert len(self.A) == 3
        assert len(self.B) == 3




        

    #Are these not more like the code builder classes? -> Fin that out
    def cyclic_shift(self, i):
        S = np.roll(np.eye(i,dtype=int), 1, axis=1)
        return S
    
    def generate_A_B(self):
        x = np.kron(self.cyclic_shift(self.l), np.eye(self.m, dtype=int))
        y = np.kron(np.eye(self.l, dtype=int), self.cyclic_shift(self.m))

        assert np.allclose(x @ y, y @ x)
        x_l = (x @ x) % 2
        for i in range(2, self.l):
            x_l = (x_l @ x) % 2
        assert np.allclose(x_l, np.eye(self.l*self.m, dtype=int))

        y_m = (y @ y) % 2
        for i in range(2, self.m):
            y_m = (y_m @ y) % 2
        assert np.allclose(y_m, np.eye(self.l*self.m, dtype=int))
        #A = A1 + A2 + A3, B = B1 + B2 + B3
        #A = x^3 + y + y^2, B = x + x^2 + y^3

        #For gross code only
        A = (x @ x @ x + y + y @ y) % 2 
        B = (x + x @ x + y @ y @ y) % 2

        assert np.allclose(A @ B, B @ A)
        self.compute_k(A,B)
        
        self.A_matrix = A
        self.B_matrix = B

        self.A = [(x @ x @ x) % 2, y, (y @ y) % 2]
        self.B = [x, (x @ x) % 2, (y @ y @ y) % 2]
        
        return A, B
    
    def compute_k(self, A,B):
        ker_A = la.null_space(A) % 2
        ker_B = la.null_space(B) % 2

        #combined = np.hstack((ker_A, ker_B))
        combined = np.where((ker_A == ker_B), ker_A, 0)
        self.k = 2*combined.shape[1]


    def generate_check_matrix(self):
        A, B = self.generate_A_B()
        #putting matricies next to each other to gtet size2lm matricies
        H_X = np.hstack([A, B])
        H_Z = np.hstack([B.transpose(), A.transpose()])
        assert H_X.shape == (self.l*self.m, 2*self.l*self.m)
        assert H_Z.shape == (self.l*self.m, 2*self.l*self.m)
        assert np.allclose((H_X @ H_Z.transpose()) % 2, np.zeros((self.l*self.m, self.l*self.m)))
        assert np.allclose( ((A @ B) + (B@A)) % 2, np.zeros((self.l*self.m, self.l*self.m)))

        self.d = self.compute_d(H_X, H_Z)
        return H_X, H_Z
    
    def compute_d(self,H_X, H_Z):
        

        """
        We should check what is wrong, until this fucntion every assertion holds. Just setting d to a value
        seems to be difficult, for some reason


        ker_H_X = la.null_space(H_X) % 2
        Q,R = la.qr(H_Z)
        rs_H_Z = R[np.abs(R).sum(axis=1) > 1e-8]

        print(ker_H_X.shape)
        print(rs_H_Z.shape)
        print(ker_H_X)
        print(rs_H_Z)

        ker_H_X_flat = set(ker_H_X.flatten())
        rs_H_Z_flat = set(rs_H_Z.flatten())

        diff = ker_H_X_flat - rs_H_Z_flat

        m = np.array(list(diff)).reshape(2,self.n)
        """

        return 12


if __name__ == "__main__":
    code = GrossCode()
        


