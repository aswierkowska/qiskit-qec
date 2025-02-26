
from qiskit_qec.codes.qec_code import QECCode
import numpy as np
import scipy.linalg as la
import sympy as sp
import belief_propagation as bp
import networkx as nx
import matplotlib.pyplot as plt
import scipy


from bposd.css import css_code
from ldpc import mod2

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

    def find_pivot_rows(self, mat):
        """
        Find the pivot rows of a given matrix.

        This function finds the pivot rows of the input matrix. The input matrix can be either a dense numpy array or 
        a sparse scipy matrix. It first converts the input matrix to a CSR format if sparse, then performs row 
        reduction to determine the pivot rows.

        Parameters
        ----------
        mat : Union[np.ndarray, scipy.sparse.spmatrix]
            The input matrix.

        Returns
        -------
        numpy.ndarray
            A 1D array of indices of the pivot rows.
        """

        # Get the shape of the matrix
        num_rows, num_cols = mat.shape

        # Track pivot rows
        pivot_rows = []
        pivot_col = 0  # Track column index for pivots

        for row in range(num_rows):
            # Find the first nonzero entry in this row at or after pivot_col
            while pivot_col < num_cols and np.all(mat[row:, pivot_col] == 0):
                pivot_col += 1  # Move to the next column

            if pivot_col >= num_cols:  # If we exceed column limits, stop
                break

            # Find the row with the largest absolute value in this column (partial pivoting)
            max_row = row + np.argmax(np.abs(mat[row:, pivot_col]))

            # Swap rows if necessary
            if max_row != row:
                mat[[row, max_row]] = mat[[max_row, row]]

            # Check if the pivot element is nonzero
            if mat[row, pivot_col] != 0:
                pivot_rows.append(row)  # Store this row as a pivot row

                # Normalize pivot row (optional, but helps stability)
                mat[row] = mat[row] / mat[row, pivot_col]

                # Eliminate all rows below the pivot
                for lower_row in range(row + 1, num_rows):
                    if mat[lower_row, pivot_col] != 0:
                        mat[lower_row] -= mat[lower_row, pivot_col] * mat[row]

            pivot_col += 1  # Move to the next column

        return np.array(pivot_rows, dtype=int)


    def gf2_rank(self,matrix):
        """
        Computes the rank of a binary matrix over GF(2) using Gaussian elimination.
        
        Parameters
        ----------
        matrix : np.ndarray
            A binary matrix (0s and 1s).

        Returns
        -------
        int
            The rank of the matrix in GF(2).
        """
        # Convert to binary (just in case the input has floating points)
        matrix = np.array(matrix, dtype=int) % 2  # Ensure everything is in GF(2)
        
        num_rows, num_cols = matrix.shape
        rank = 0
        pivot_col = 0  # Tracks column index for pivots

        for row in range(num_rows):
            # Find a pivot column (first nonzero entry)
            while pivot_col < num_cols and np.all(matrix[row:, pivot_col] == 0):
                pivot_col += 1  # Move to the next column

            if pivot_col >= num_cols:  # No more pivots possible
                break

            # Find the row with a 1 in the pivot_col and swap
            max_row = row + np.argmax(matrix[row:, pivot_col])  # Finds first 1 in column
            matrix[[row, max_row]] = matrix[[max_row, row]]  # Swap rows

            # Zero out all other 1s in this column
            for lower_row in range(row + 1, num_rows):
                if matrix[lower_row, pivot_col] == 1:
                    matrix[lower_row] ^= matrix[row]  # XOR the row (mod 2 subtraction)

            rank += 1
            pivot_col += 1  # Move to the next column

        return rank



    def row_echelon(self,matrix, full=False):

        num_rows, num_cols = np.shape(matrix)

        # Take copy of matrix if numpy (why?) and initialise transform matrix to identity
        if isinstance(matrix, np.ndarray):
            the_matrix = np.copy(matrix)
            transform_matrix = np.identity(num_rows).astype(int)
        elif isinstance(matrix, scipy.sparse.csr.csr_matrix):
            the_matrix = matrix
            transform_matrix = scipy.sparse.eye(num_rows, dtype="int", format="csr")
        else:
            raise ValueError("Unrecognised matrix type")

        pivot_row = 0
        pivot_cols = []

        # Iterate over cols, for each col find a pivot (if it exists)
        for col in range(num_cols):
            # Select the pivot - if not in this row, swap rows to bring a 1 to this row, if possible
            if the_matrix[pivot_row, col] != 1:
                # Find a row with a 1 in this col
                swap_row_index = pivot_row + np.argmax(the_matrix[pivot_row:num_rows, col])

                # If an appropriate row is found, swap it with the pivot. Otherwise, all zeroes - will loop to next col
                if the_matrix[swap_row_index, col] == 1:
                    # Swap rows
                    the_matrix[[swap_row_index, pivot_row]] = the_matrix[
                        [pivot_row, swap_row_index]
                    ]

                    # Transformation matrix update to reflect this row swap
                    transform_matrix[[swap_row_index, pivot_row]] = transform_matrix[
                        [pivot_row, swap_row_index]
                    ]

            # If we have got a pivot, now let's ensure values below that pivot are zeros
            if the_matrix[pivot_row, col]:
                if not full:
                    elimination_range = [k for k in range(pivot_row + 1, num_rows)]
                else:
                    elimination_range = [k for k in range(num_rows) if k != pivot_row]

                # Let's zero those values below the pivot by adding our current row to their row
                for j in elimination_range:
                    if (
                        the_matrix[j, col] != 0 and pivot_row != j
                    ):  ### Do we need second condition?
                        the_matrix[j] = (the_matrix[j] + the_matrix[pivot_row]) % 2

                        # Update transformation matrix to reflect this op
                        transform_matrix[j] = (
                            transform_matrix[j] + transform_matrix[pivot_row]
                        ) % 2

                pivot_row += 1
                pivot_cols.append(col)

            # Exit loop once there are no more rows to search
            if pivot_row >= num_rows:
                break

        # The rank is equal to the maximum pivot index
        matrix_rank = pivot_row
        row_esch_matrix = the_matrix

        return [row_esch_matrix, matrix_rank, transform_matrix, pivot_cols]



    def gf2_nullspace(self,matrix):

        transpose = matrix.T
        m, n = transpose.shape
        _, matrix_rank, transform, _ = self.row_echelon(transpose)
        nspace = transform[matrix_rank:m]
        return nspace


    def gf2_pivot_rows(self,matrix):
        """
        Finds the pivot rows of a binary matrix over GF(2) using Gaussian elimination.

        Parameters
        ----------
        matrix : np.ndarray
            A binary matrix (0s and 1s).

        Returns
        -------
        np.ndarray
            An array containing the indices of the pivot rows.
        """
        # Convert to binary (mod 2) to ensure correctness
        matrix = np.array(matrix, dtype=int) % 2
        num_rows, num_cols = matrix.shape

        pivot_rows = []
        pivot_col = 0  # Track the pivot column index

        for row in range(num_rows):
            # Find a pivot column (first nonzero entry in a column)
            while pivot_col < num_cols and np.all(matrix[row:, pivot_col] == 0):
                pivot_col += 1  # Move to the next column

            if pivot_col >= num_cols:  # No more pivots possible
                break

            # Find the first row with a 1 in pivot_col
            max_row = row + np.argmax(matrix[row:, pivot_col])  # Finds first 1 in column
            matrix[[row, max_row]] = matrix[[max_row, row]]  # Swap rows

            # Store the pivot row index
            pivot_rows.append(row)

            # Zero out all rows below the pivot
            for lower_row in range(row + 1, num_rows):
                if matrix[lower_row, pivot_col] == 1:
                    matrix[lower_row] ^= matrix[row]  # XOR the row (mod 2)

            pivot_col += 1  # Move to the next column

        return np.array(pivot_rows, dtype=int)  # Return pivot row indices





    def compute_logicals(self,hx,hz):
        hx = scipy.sparse.csr_matrix(hx)
        hz = scipy.sparse.csr_matrix(hz)
        ker_hx = scipy.sparse.csr_matrix(self.gf2_nullspace(hx.toarray()))
        log_stack = scipy.sparse.vstack([hz, ker_hx])

        rank_hz = self.gf2_rank(hz.toarray())
        pivots = mod2.pivot_rows(log_stack)[rank_hz:]
        #pirvot_v2 = self.find_pivot_rows(log_stack.toarray())[rank_hz:]
        #assert np.allclose(pivots, pirvot_v2)
        log_ops = log_stack[pivots]

        rank_hz_v2 = mod2.rank(hz)
        assert np.allclose(rank_hz, rank_hz_v2)
        ker_hx_v2 = self.gf2_nullspace(hx.toarray())
        assert np.allclose(ker_hx.toarray(), ker_hx_v2)
        pivots_v2 = self.gf2_pivot_rows(log_stack.toarray())[rank_hz:]
        print("Pivots: ", pivots)
        print("Pivots V2: ", pivots_v2)
        assert np.allclose(pivots, pivots_v2)

        return log_ops

    

    def _logical_z(self, n):
        logz = self.null_space_mod2(self.H_X)
        z_log = []
        for i in logz:
            z_log.append(self.indicator_vector(i,n))
        

        z_log = self.compute_logicals(self.H_X, self.H_Z)
        z_log = z_log.toarray()
        return z_log

    def _logical_x(self, n):
        logx = self.null_space_mod2(self.H_Z)
        x_log = []
        for i in logx:
            x_log.append(self.indicator_vector(i,n))

        x_log = self.compute_logicals(self.H_Z, self.H_X)
        x_log = x_log.toarray()
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

    qcode = css_code(code.H_X, code.H_Z)
    lx = qcode.lx
    lz = qcode.lz

    print("Logical X: ", lx)
    print(" Code Logical X: ", code.logical_x)


