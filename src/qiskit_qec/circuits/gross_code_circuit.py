"""Generates circuits for the gross code."""
from qiskit_qec.circuits.code_circuit import CodeCircuit
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

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

        self.dept_8_syndrome_measurement()
    

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


    #prpoably we need these functions
    #all thesApply a logical Z gate to the code.e are listed in the syndrome measurenemnt
    def CNOT(self, c, t):
        """CNOT with control qubit c and taget qubit t"""
        pass

    def InitX(self):
        """Initialize qubit q in the state |+> = (|0> + |1>)/sqrt(2)"""
        pass

    def InitZ(self,q):
        """Initialize qubit q in the state |0>"""
        pass

    def MeasX(self,q):
        """Measure qubit q in the X basis, |+> or |->"""
        pass

    def MeasZ(self,q):
        """Measure qubit q in the Z basis, |0> or |1>"""
        pass

    def Idle(self,q):
        """Idle operation, Identity on qubit q"""
        pass

    #helper
    def _get_j(A, i):
        """
        Get j for the given i
        """
        #return column with the 1 in row i
        return list(A[i]).index(1)

    def dept_8_syndrome_measurement(self):
        """
        Depth-8 syndrom measurenment cycle circuit
        """
        j = 0 #figure out what j is and where it comes from
        #Round 1
        for i in range(self.n/2):
            self.InitX(self.qr_X[i])
            self.CNOT(self.qr_right[self._get_j(self.A[0].transpose(), i)], self.qr_Z[i])
            self.Idle(self.qr_left[i])

        #Round 2
        for i in range(self.n/2):
            self.CNOT(self.qr_X[i], self.qr_left[self._get_j(self.A[1],i )])
            self.CNOT(self.qr_right[self._get_j(self.A[2].transpose(),i)], self.qr_Z[i])

        #Round 3
        for i in range(self.n/2):
            self.CNOT(self.qr_X[i], self.qr_right[self._get_j(self.B[2],i)])
            self.CNOT(self.qr_left[self._get_j(self.B[0].transpose(), i)], self.qr_Z[i])

        #Round 4
        for i in range(self.n/2):
            self.CNOT(self.qr_X[i], self.qr_right[self._get_j(self.B[0], i)])
            self.CNOT(self.qr_left[self._get_j(self.B[2].transpose(), i)], self.qr_Z[i])

        #Round 5
        for i in range(self.n/2):
            self.CNOT(self.qr_X[i], self.qr_right[self._get_j(self.B[2], i)])
            self.CNOT(self.qr_left[self._get_j(self.B[2].transpose(), i)], self.qr_Z[i])

        #Round 6
        for i in range(self.n/2):
            self.CNOT(self.qr_X[i], self.qr_left[self._get_j(self.A[0], i)])
            self.CNOT(self.qr_right[self._get_j(self.A[1].transpose(), i)], self.qr_Z[i])

        #Round 7
        for i in range(self.n/2):
            self.CNOT(self.qr_X[i], self.qr_left[self._get_j(self.A[2], i)])
            self.MeasZ(self.qr_Z[i])
            self.Idle(self.qr_right[i])

        #Round 8
        for i in range(self.n/2):
            self.MeasX(self.qr_X[i])
            self.InitZ(self.qr_Z[i])
            self.Idle(self.qr_left[i])
            self.Idle(self.qr_right[i])


    def readout(self):
        """
        Readout of all logical qubits
        """
        pass

if __name__ == "__main__":
    code = GrossCode()
    circuit = GrossCodeCircuit(code)
    pass