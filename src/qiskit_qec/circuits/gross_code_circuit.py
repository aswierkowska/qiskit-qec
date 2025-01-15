"""Generates circuits for the gross code."""
from qiskit_qec.circuits.code_circuit import CodeCircuit
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
import sys
print(sys.path)


from qiskit_qec.codes.gross_code import GrossCode

class GrossCodeCircuit(CodeCircuit):
    
    #TODO: Figure out what inputs are needed
    def __init__(self, code, T: int = 7):
        super().__init__()
        self.code = code
        self._get_code_properties()
        self.T = T

        self.circuit = None
        states = ["0", "1"]

        #add one point we need that
        self.circuit = QuantumCircuit()

        """The SM circuit uses 2n physical qubits in total: n data qubits and n ancillary check qubits. Page 12 Quantum memeory paper"""
        assert self.n % 2 == 0
        self.qr_left = QuantumRegister(int(self.n/2), 'qr_left')
        self.qr_right = QuantumRegister(int(self.n/2), 'qr_right')
        self.qr_X = QuantumRegister(int(self.n/2), 'qr_X')
        self.qr_Z = QuantumRegister(int(self.n/2), 'qr_Z')

        self.cr_X = ClassicalRegister(int(self.n/2), 'cr_X')
        self.cr_Z = ClassicalRegister(int(self.n/2), 'cr_Z')

        self.circuit.add_register(self.qr_X)
        self.circuit.add_register(self.qr_left)
        self.circuit.add_register(self.qr_right)
        self.circuit.add_register(self.qr_Z)

        self.circuit.add_register(self.cr_X)
        self.circuit.add_register(self.cr_Z)

        self.connectivity_dict_L = {}
        self.connectivity_dict_R = {}
        self._init_connectivity_dict()

    

    def _get_code_properties(self):
        """ Can this be done nicer?"""
        self.n = self.code.n
        self.k = self.code.k
        self.d = self.code.d
        self.A = self.code.A
        self.B = self.code.B
        self.H_X = self.code.H_X
        self.H_Z = self.code.H_Z
        self.A_matrix = self.code.A_matrix
        self.B_matrix = self.code.B_matrix
    
    def _init_connectivity_dict(self):
        """Initialize the connectivity dict"""
        for i in range(int(self.n/2)):
            self.connectivity_dict_L[i] = set([])
            self.connectivity_dict_R[i] = set([])


    #prpoably we need these functions
    #all thesApply a logical Z gate to the code.e are listed in the syndrome measurenemnt
    def CNOT(self, c, t):
        """CNOT with control qubit c and taget qubit t"""
        self.circuit.cx(c, t)

    def InitX(self,q):
        """Initialize qubit q in the state |+> = (|0> + |1>)/sqrt(2)"""
        self.circuit.initialize([1,0],q)
        self.circuit.h(q)


    def InitZ(self,q):
        """Initialize qubit q in the state |0>
        Qubits are already initiliazed to |0> in qiskit
        But mabye for second round we need to do this?
        """
        self.circuit.initialize([1,0],q)
        pass

    def MeasX(self,q, c):
        """Measure qubit q in the X basis, |+> or |->"""
        self.circuit.h(q)
        self.circuit.measure(q, c) # Figure out how to do this

    def MeasZ(self,q, c):
        """Measure qubit q in the Z basis, |0> or |1>"""
        self.circuit.measure(q, c) # Figure out how to do this
        pass

    def Idle(self,q):
        """Idle operation, Identity on qubit q"""
        pass

    #helper
    def _get_j(self,A, i):
        """
        Get j for the given i
        """
        #return column with the 1 in row i
        return list(A[i]).index(1)
    
    def add_barrier(self):
        """Add a barrier to the circuit"""
        self.circuit.barrier()

    def depth_8_syndrome_measurement(self):
        """
        Depth-8 syndrom measurenment cycle circuit
        """
        j = 0 #figure out what j is and where it comes from
        L = 8
        R = 27
        #Round 1
        for i in range(int(self.n/2)):
            if self._get_j(self.A[0].transpose(), i) == R: print("R connects to Z: ", i)

            self.InitX(self.qr_X[i])
            self.CNOT(self.qr_right[self._get_j(self.A[0].transpose(), i)], self.qr_Z[i])
            self.Idle(self.qr_left[i])

            self.connectivity_dict_R[i].add(self._get_j(self.A[0].transpose(), i))
        
        print("Round 1 done")
        self.add_barrier()

        #Round 2
        for i in range(int(self.n/2)):
            if self._get_j(self.A[1], i) == L: print("L connects to X: ", i)
            if self._get_j(self.A[2].transpose(), i) == R: print("R connects to Z: ", i)
                
            self.CNOT(self.qr_X[i], self.qr_left[self._get_j(self.A[1],i )])
            self.CNOT(self.qr_right[self._get_j(self.A[2].transpose(),i)], self.qr_Z[i])

            self.connectivity_dict_L[i].add(self._get_j(self.A[1], i))
            self.connectivity_dict_R[i].add(self._get_j(self.A[2].transpose(), i))
        
        print("Round 2 done")
        self.add_barrier()

        #Round 3
        for i in range(int(self.n/2)):
            if self._get_j(self.B[1], i) == R: print("R connects to X: ", i)
            if self._get_j(self.B[0].transpose(), i) == L: print("L connects to Z: ", i)

            self.CNOT(self.qr_X[i], self.qr_right[self._get_j(self.B[1],i)])
            self.CNOT(self.qr_left[self._get_j(self.B[0].transpose(), i)], self.qr_Z[i])

            self.connectivity_dict_R[i].add(self._get_j(self.B[1], i))
            self.connectivity_dict_L[i].add(self._get_j(self.B[0].transpose(), i))
        
        print("Round 3 done")
        self.add_barrier()

        #Round 4
        for i in range(int(self.n/2)):
            if self._get_j(self.B[0], i) == R: print("R connects to X: ", i)
            if self._get_j(self.B[1].transpose(), i) == L: print("L connects to Z: ", i)

            self.CNOT(self.qr_X[i], self.qr_right[self._get_j(self.B[0], i)])
            self.CNOT(self.qr_left[self._get_j(self.B[1].transpose(), i)], self.qr_Z[i])

            self.connectivity_dict_R[i].add(self._get_j(self.B[0], i))
            self.connectivity_dict_L[i].add(self._get_j(self.B[1].transpose(), i))
        
        print("Round 4 done")
        self.add_barrier()

        #Round 5
        for i in range(int(self.n/2)):
            if self._get_j(self.B[2], i) == R: print("R connects to X: ", i)
            if self._get_j(self.B[2].transpose(), i) == L: print("L connects to Z: ", i)

            self.CNOT(self.qr_X[i], self.qr_right[self._get_j(self.B[2], i)])
            self.CNOT(self.qr_left[self._get_j(self.B[2].transpose(), i)], self.qr_Z[i])

            self.connectivity_dict_R[i].add(self._get_j(self.B[2], i))
            self.connectivity_dict_L[i].add(self._get_j(self.B[2].transpose(), i))
        
        print("Round 5 done")
        self.add_barrier()

        #Round 6
        for i in range(int(self.n/2)):
            if self._get_j(self.A[0], i) == L: print("L connects to X: ", i)
            if self._get_j(self.A[1].transpose(), i) == R: print("R connects to Z: ", i)

            self.CNOT(self.qr_X[i], self.qr_left[self._get_j(self.A[0], i)])
            self.CNOT(self.qr_right[self._get_j(self.A[1].transpose(), i)], self.qr_Z[i])

            self.connectivity_dict_L[i].add(self._get_j(self.A[0], i))
            self.connectivity_dict_R[i].add(self._get_j(self.A[1].transpose(), i))
        
        print("Round 6 done")
        self.add_barrier()

        #Round 7
        for i in range(int(self.n/2)):
            if self._get_j(self.A[2], i) == L: print("L connects to X: ", i)

            self.CNOT(self.qr_X[i], self.qr_left[self._get_j(self.A[2], i)])
            self.MeasZ(self.qr_Z[i], self.cr_Z[i])
            self.Idle(self.qr_right[i])

            self.connectivity_dict_L[i].add(self._get_j(self.A[2], i))

        print("Round 7 done")
        self.add_barrier()

        #Round 8
        for i in range(int(self.n/2)):
            self.MeasX(self.qr_X[i], self.cr_X[i])
            self.InitZ(self.qr_Z[i])
            self.Idle(self.qr_left[i])
            self.Idle(self.qr_right[i])

        self.add_barrier()

    
    def multiple_syndrome_measurement(self, T):
        """
        Multiple syndrome measurement cycle circuit
        """
        for i in range(T):
            self.circuit.barrier()
            self.depth_8_syndrome_measurement()


    def readout(self):
        """
        Readout of all logical qubits
        """
        pass

    def verify_connectivity(self):
        """
        Verify that the connectivity is correct
        """
        #check length for each value
        for k,v in enumerate(self.connectivity_dict_L):
            assert len(self.connectivity_dict_L[k]) == 6
        
        for k,v in enumerate(self.connectivity_dict_R):
            assert len(self.connectivity_dict_R[k]) == 6

        print("Connectivity is correct")

    #We may need to check out this fuction again becuase i have no idea what it does
    def check_nodes(self, nodes, ignore_extras=False, minimal=False):
        return super().check_nodes(nodes, ignore_extras, minimal)
    
    def is_cluster_neutral(self, nodes):
        return super().is_cluster_neutral(nodes)
    
    def string2nodes(self, string, **kwargs):
        return super().string2nodes(string, **kwargs)

if __name__ == "__main__":
    code = GrossCode()
    circuit = GrossCodeCircuit(code)
    circuit.multiple_syndrome_measurement(2)
    #circuit.circuit.draw(output='mpl', filename='gross_code_circuit.png', vertical_compression='high', scale=0.3, fold=500)
    circuit.verify_connectivity()
    pass