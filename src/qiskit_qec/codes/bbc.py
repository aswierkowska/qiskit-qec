
from qiskit_qec.codes.qec_code import QECCode
import numpy as np
import scipy.linalg as la
import sympy as sp
import belief_propagation as bp
import networkx as nx
import matplotlib.pyplot as plt

class BBCode(QECCode):

    def __init__(self, n=144, k=12, distance=12,l=12,m=6,exp_A=[3,1,2], exp_B=[3,1,2]):
        """
        Args:
            n: The code length
            k: The number of logical qubits in the code.
            d: The distance of the code.
            exp_A: Exponents of the A matrix
            exp_B: Exponents of the B matrix
        """
        assert len(exp_A) == 3
        assert len(exp_B) == 3
        super().__init__()


        self.l = l
        self.m = m

        self.n = 2*self.l*self.m
        self.k = k
        self.d = distance

        self.A = []
        self.B = []

        self.A_matrix = None
        self.B_matrix = None

        self.H_X, self.H_Z = self.generate_check_matrix(exp_A,exp_B)
        self.H = np.vstack([self.H_X, self.H_Z])
        #parameters for gross code
        #assert self.k == k
        assert self.d == distance
        assert self.n == n
        assert self.H_X.shape == (self.l * self.m, self.n)
        assert self.H_Z.shape == (self.l * self.m, self.n)
        assert len(self.A) == 3
        assert len(self.B) == 3


        self.x_gauges = self._x_gauge_geneators(self.d)
        self.z_gauges = self._z_gauge_geneators(self.d)
        self.x_stabilizers = self._x_stabilizer_geneators(self.d)
        self.z_stabilizers = self._z_stabilizer_geneators(self.d)
        self.logical_z = self._logical_z(self.n)
        self.logical_x = self._logical_x(self.n)

        #print(self.logical_x)
        #print(self.logical_z)



    def get_indices_of_ones(self,row):
        """
        Returns the indices of all entries in the row that are equal to 1.
        """
        return [index for index, value in enumerate(row) if value == 1]
    

    def indicator_vector(self,positions, n):
        """Creates a binary vector of length n with 1s at the specified positions."""
        vector = np.zeros(n, dtype=int)  # Initialize all zeros
        vector[positions] = 1  # Set the given indices to 1
        return vector

        

    #TODO thses functions but how? should not make sense i guess?
    #Because gross code is more like a surface code this should be fine. It's not CSS Code tho
    def _x_gauge_geneators(self, d):
        x_gauges = []
        for i in range(int(self.n/2)):
            l = self.get_indices_of_ones(self.H_X[i])
            x_gauges.append([l[0],l[3],l[4],l[5],l[1],l[2]])
            #x_gauges.append(self.get_indices_of_ones(self.H_X[i]))
        return x_gauges        

    def _z_gauge_geneators(self, d):
        z_gauges = []
        for i in range(int(self.n/2)):
            l = self.get_indices_of_ones(self.H_Z[i])
            z_gauges.append([l[0],l[3],l[4],l[5],l[1],l[2]])
            #z_gauges.append(self.get_indices_of_ones(self.H_Z[i]))
        return z_gauges
    
    def _x_stabilizer_geneators(self, d):
        x_stabilizers = []
        for i in range(self.H_X.shape[0]):
            x_stabilizers.append(self.get_indices_of_ones(self.H_X[i]))
        return x_stabilizers


    def _z_stabilizer_geneators(self, d):
        z_stabilizers = []
        for i in range(self.H_Z.shape[0]):
            z_stabilizers.append(self.get_indices_of_ones(self.H_Z[i]))
        return z_stabilizers
    
    def null_space_mod2(self, matrix):
        M = sp.Matrix(matrix)
        M = M.rref(iszerofunc=lambda x: x % 2 == 0)[0]

        null_vectors = M.nullspace()
        null_vectors_mod2 = [np.array(v) % 2 for v in null_vectors]
        logicals = [list(np.where(v == 1)[0]) for v in null_vectors_mod2]
        return logicals

    def _logical_z(self, n):
        logz = self.null_space_mod2(self.H_X)
        z_log = []
        for i in logz:
            z_log.append(self.indicator_vector(i,n))
        return z_log

    def _logical_x(self, n):
        logx = self.null_space_mod2(self.H_Z)
        x_log = []
        for i in logx:
            x_log.append(self.indicator_vector(i,n))
        return x_log

    def cyclic_shift(self, i):
        S = np.roll(np.eye(i,dtype=int), 1, axis=1)
        return S
    

    def matrix_exponentiation(self, m, exp):
        assert exp >= 0
        result = np.eye(m.shape[0], dtype=int)  # Identity matrix
        base = m  # Copy of the matrix

        while exp > 0:
            if exp % 2 == 1:  # If exp is odd, multiply result by base
                result = result @ base % 2
            base = base @ base % 2 # Square the base
            exp //= 2  # Reduce exponent by half

        return result % 2  # Apply modulo 2

    
    def generate_A_B(self,exp_A, exp_B):
        x = np.kron(self.cyclic_shift(self.l), np.eye(self.m, dtype=int))
        y = np.kron(np.eye(self.l, dtype=int), self.cyclic_shift(self.m))

        assert np.allclose((x @ y) % 2, (y @ x) % 2)
        x_l = (x @ x) % 2
        for i in range(2, self.l):
            x_l = (x_l @ x) % 2
        assert np.allclose(x_l, np.eye(self.l*self.m, dtype=int))

        y_m = (y @ y) % 2
        for i in range(2, self.m):
            y_m = (y_m @ y) % 2
        assert np.allclose(y_m, np.eye(self.l*self.m, dtype=int))
        assert np.allclose(x_l, y_m)
        #A = A1 + A2 + A3, B = B1 + B2 + B3
        #A = x^3 + y + y^2, B = x + x^2 + y^3

        #For gross code only
        #A = (x @ x @ x + y + y @ y) % 2 
        #B = (x + x @ x + y @ y @ y) % 2

        A = (self.matrix_exponentiation(x, exp_A[0]) + self.matrix_exponentiation(y, exp_A[1]) + self.matrix_exponentiation(y, exp_A[2])) % 2
        B = (self.matrix_exponentiation(y, exp_B[0]) + self.matrix_exponentiation(x, exp_B[1]) + self.matrix_exponentiation(x, exp_B[2])) % 2

        assert np.allclose(A @ B, B @ A)
        self.check_A_B_properties(A,B)
        #self.compute_k(A,B)
        
        self.A_matrix = A
        self.B_matrix = B

        #self.A = [(x @ x @ x) % 2, y, (y @ y) % 2]
        #self.B = [x, (x @ x) % 2, (y @ y @ y) % 2]
        self.A = [self.matrix_exponentiation(x, exp_A[0]), self.matrix_exponentiation(y, exp_A[1]), self.matrix_exponentiation(y, exp_A[2])]
        self.B = [self.matrix_exponentiation(y, exp_B[0]), self.matrix_exponentiation(x, exp_B[1]), self.matrix_exponentiation(x, exp_B[2])]
        

        return A, B
    
    def check_A_B_properties(self,A,B):
        #each row sums up to 3
        #each column sums up to 3
        assert np.allclose(A.sum(axis=1), 3*np.ones(A.shape[0]))
        assert np.allclose(A.sum(axis=0), 3*np.ones(A.shape[1]))
        assert np.allclose(B.sum(axis=1), 3*np.ones(B.shape[0]))
        assert np.allclose(B.sum(axis=0), 3*np.ones(B.shape[1]))
        print("A and B properties are satisfied")
    
    def compute_k(self, A,B):
        ker_A = la.null_space(A) % 2
        ker_B = la.null_space(B) % 2

        combined = np.where((ker_A == ker_B), ker_A, 0)
        self.k = 2*combined.shape[1]


    def generate_check_matrix(self, exp_A, exp_B):
        A, B = self.generate_A_B(exp_A, exp_B)
        #putting matricies next to each other to gtet size2lm matricies
        H_X = np.hstack([A, B])
        H_Z = np.hstack([B.transpose(), A.transpose()])
        assert H_X.shape == (self.l*self.m, 2*self.l*self.m)
        assert H_Z.shape == (self.l*self.m, 2*self.l*self.m)
        assert np.allclose((H_X @ H_Z.transpose()) % 2, np.zeros((self.l*self.m, self.l*self.m)))
        assert np.allclose( ((A @ B) + (B@A)) % 2, np.zeros((self.l*self.m, self.l*self.m)))
        #each row sums up to 6
        #each column sums up to 6
        assert np.allclose(H_X.sum(axis=1), 6*np.ones(H_X.shape[0]))
        assert np.allclose(H_Z.sum(axis=1), 6*np.ones(H_Z.shape[0]))

        self.d = self.compute_d(H_X, H_Z)
        k = 2 * self.l * self.m - np.linalg.matrix_rank(H_X) - np.linalg.matrix_rank(H_Z)
        return H_X, H_Z
    
    def compute_d(self,H_X, H_Z):
        

        """
        We should check what is wrong, until this fucntion every assertion holds. Just setting d to a value
        seems to be difficult, for some reason


        ker_H_X = la.null_space(H_X) % 2
        Q,R = la.qr(H_Z)
        rs_H_Z = R[np.abs(R).sum(axis=1) > 1e-8]

        ker_H_X_flat = set(ker_H_X.flatten())
        rs_H_Z_flat = set(rs_H_Z.flatten())

        diff = ker_H_X_flat - rs_H_Z_flat

        m = np.array(list(diff)).reshape(2,self.n)
        """

        return self.d
    
    def verify_syndrome(self):
        H = self.H
        error = np.zeros(self.n, dtype=int)
        error[0] = 1
        error[1] = 1
        syndrome = (error @ H) % 2
        print("Syndrom: ", syndrome)
        print("Error: ", error)
        print(H[0])
        print(H[1])

        pass

    def show_tannergraph(self):
        tg = bp.TannerGraph.from_biadjacency_matrix(self.H,channel_model=1)
        g = tg.to_nx()

        # Assuming `g` is your Tanner graph
        fig = plt.figure()

        # Get the two bipartite sets
        top_set, bottom_set = nx.bipartite.sets(g)   # Remaining nodes are in `bottom`

        # Sort nodes to ensure consistent ordering
        top = sorted(top_set)
        bottom = sorted(bottom_set)

        # Assign colors
        node_colors = {}

        # Color first 72 check nodes differently from the next 72
        for i, node in enumerate(top):
            node_colors[node] = "red" if i < int(self.n/2) else "blue"

        # Color first 72 variable nodes differently from the next 72
        for i, node in enumerate(bottom):
            node_colors[node] = "green" if i < int(self.n/2) else "orange"

        # Get color list for drawing
        node_color_list = [node_colors[n] for n in g.nodes()]

        # Labels for nodes
        labels = {node: d["label"] for node, d in g.nodes(data=True)}

        # Draw the graph
        nx.draw(
            g,
            pos=nx.bipartite_layout(g, top),  # Layout for bipartite graph
            with_labels=True,
            labels=labels,
            node_color=node_color_list
        )

        plt.show()


if __name__ == "__main__":
    code = BBCode(72,12,6,6,6,[3,1,2],[3,1,2])
    #np.set_printoptions(threshold=np.inf)
    #code._logical_x(144)
    #code.show_tannergraph()
    code.verify_syndrome()


