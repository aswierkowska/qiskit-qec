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

        self.qc = None
        self.circuit = {}
        self.noisy_circuit = {}
        states = ["0", "1"]
        self.basis = "z"

        #add one point we need that
        self.qc = QuantumCircuit()

        """The SM circuit uses 2n physical qubits in total: n data qubits and n ancillary check qubits. Page 12 Quantum memeory paper"""
        assert self.n % 2 == 0
        self.n_half = int(self.n/2)#just for convinience
        self.qr_X = QuantumRegister(int(self.n/2), 'qr_X')
        self.qr_Z = QuantumRegister(int(self.n/2), 'qr_Z')
        self.qr_left_right = QuantumRegister(int(self.n), 'qr_left_right')
        self.cr_final = ClassicalRegister(int(self.n), 'final_readout')

        #Register Layout taken from the paper
        self.qc.add_register(self.qr_X)
        self.qc.add_register(self.qr_left_right)
        self.qc.add_register(self.qr_Z)
        self.qc.add_register(self.cr_final)

        self.connectivity_dict_L = {}
        self.connectivity_dict_R = {}
        self._init_connectivity_dict()


        self._gauges4stabilizers = []
        self._stabilizers = [self.x_stabilizers, self.z_stabilizers]
        self._gauges = [self.x_gauges, self.z_gauges]
        for j in range(2):
            self._gauges4stabilizers.append([])
            for stabilizer in self._stabilizers[j]:
                gauges = []
                for g, gauge in enumerate(self._gauges[j]):
                    if set(stabilizer).intersection(set(gauge)) == set(gauge):
                        gauges.append(g)
                self._gauges4stabilizers[j].append(gauges)

        self.multiple_syndrome_measurement(self.T)

        for state in states:
            self.circuit[state] = self.qc.copy()
        
        self.detectors, self.logicals = self.stim_detectors()

    

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
        self.x_gauges = self.code.x_gauges
        self.z_gauges = self.code.z_gauges
        self.x_stabilizers = self.code.x_stabilizers
        self.z_stabilizers = self.code.z_stabilizers
        self.logical_x = self.code.logical_x
        self.logical_z = self.code.logical_z


    
    def _init_connectivity_dict(self):
        """Initialize the connectivity dict"""
        for i in range(int(self.n/2)):
            self.connectivity_dict_L[i] = set([])
            self.connectivity_dict_R[i] = set([])


    #prpoably we need these functions
    #all thesApply a logical Z gate to the code.e are listed in the syndrome measurenemnt
    def CNOT(self, c, t):
        """CNOT with control qubit c and taget qubit t"""
        self.qc.cx(c, t)

    def InitX(self,q):
        """Initialize qubit q in the state |+> = (|0> + |1>)/sqrt(2)"""
        self.qc.initialize('0',q)
        self.qc.h(q)


    def InitZ(self,q):
        """Initialize qubit q in the state |0>
        Qubits are already initiliazed to |0> in qiskit
        But mabye for second round we need to do this?
        """
        self.qc.initialize('0',q)
        pass

    def MeasX(self,q, c):
        """Measure qubit q in the X basis, |+> or |->"""
        self.qc.h(q)
        self.qc.measure(q, c) # Figure out how to do this
        self.qc.reset(q)

    def MeasZ(self,q, c):
        """Measure qubit q in the Z basis, |0> or |1>"""
        
        self.qc.measure(q, c) # Figure out how to do this
        self.qc.reset(q)
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
        self.qc.barrier()

    def depth_8_syndrome_measurement(self,t):
        """
        Depth-8 syndrom measurenment cycle circuit
        """
        reg_Z = f"round_{t}_z_bits"
        cr_z = ClassicalRegister(int(self.n/2), name=reg_Z)
        self.qc.add_register(cr_z)

        reg_X = f"round_{t}_x_bits"
        #self.qc.h(q) is this necessary?
        cr_x = ClassicalRegister(int(self.n/2), name=reg_X)
        self.qc.add_register(cr_x)

        #Round 1
        for i in range(int(self.n/2)):
            #if self._get_j(self.A[0].transpose(), i) == R: print("R connects to Z: ", i)

            self.InitX(self.qr_X[i])
            self.CNOT(self.qr_left_right[self._get_j(self.A[0].transpose(), i) + self.n_half], self.qr_Z[i])
            self.Idle(self.qr_left_right[i])

            self.connectivity_dict_R[i].add(self._get_j(self.A[0].transpose(), i))
        
        self.add_barrier()

        #Round 2
        for i in range(int(self.n/2)):
            #if self._get_j(self.A[1], i) == L: print("L connects to X: ", i)
            #if self._get_j(self.A[2].transpose(), i) == R: print("R connects to Z: ", i)
                
            self.CNOT(self.qr_X[i], self.qr_left_right[self._get_j(self.A[1],i )])
            self.CNOT(self.qr_left_right[self._get_j(self.A[2].transpose(),i) + self.n_half], self.qr_Z[i])

            self.connectivity_dict_L[i].add(self._get_j(self.A[1], i))
            self.connectivity_dict_R[i].add(self._get_j(self.A[2].transpose(), i))
        
        self.add_barrier()

        #Round 3
        for i in range(int(self.n/2)):
            #if self._get_j(self.B[1], i) == R: print("R connects to X: ", i)
            #if self._get_j(self.B[0].transpose(), i) == L: print("L connects to Z: ", i)

            self.CNOT(self.qr_X[i], self.qr_left_right[self._get_j(self.B[1],i) + self.n_half])
            self.CNOT(self.qr_left_right[self._get_j(self.B[0].transpose(), i)], self.qr_Z[i])

            self.connectivity_dict_R[i].add(self._get_j(self.B[1], i))
            self.connectivity_dict_L[i].add(self._get_j(self.B[0].transpose(), i))
        
        self.add_barrier()

        #Round 4
        for i in range(int(self.n/2)):
            #if self._get_j(self.B[0], i) == R: print("R connects to X: ", i)
            #if self._get_j(self.B[1].transpose(), i) == L: print("L connects to Z: ", i)

            self.CNOT(self.qr_X[i], self.qr_left_right[self._get_j(self.B[0], i) + self.n_half])
            self.CNOT(self.qr_left_right[self._get_j(self.B[1].transpose(), i)], self.qr_Z[i])

            self.connectivity_dict_R[i].add(self._get_j(self.B[0], i))
            self.connectivity_dict_L[i].add(self._get_j(self.B[1].transpose(), i))
        
        self.add_barrier()

        #Round 5
        for i in range(int(self.n/2)):
            #if self._get_j(self.B[2], i) == R: print("R connects to X: ", i)
            #if self._get_j(self.B[2].transpose(), i) == L: print("L connects to Z: ", i)

            self.CNOT(self.qr_X[i], self.qr_left_right[self._get_j(self.B[2], i) + self.n_half])
            self.CNOT(self.qr_left_right[self._get_j(self.B[2].transpose(), i)], self.qr_Z[i])

            self.connectivity_dict_R[i].add(self._get_j(self.B[2], i))
            self.connectivity_dict_L[i].add(self._get_j(self.B[2].transpose(), i))
        
        self.add_barrier()

        #Round 6
        for i in range(int(self.n/2)):
            #if self._get_j(self.A[0], i) == L: print("L connects to X: ", i)
            #if self._get_j(self.A[1].transpose(), i) == R: print("R connects to Z: ", i)

            self.CNOT(self.qr_X[i], self.qr_left_right[self._get_j(self.A[0], i)])
            self.CNOT(self.qr_left_right[self._get_j(self.A[1].transpose(), i) + self.n_half], self.qr_Z[i])

            self.connectivity_dict_L[i].add(self._get_j(self.A[0], i))
            self.connectivity_dict_R[i].add(self._get_j(self.A[1].transpose(), i))
        
        self.add_barrier()

        #Round 7
        for i in range(int(self.n/2)):
            #if self._get_j(self.A[2], i) == L: print("L connects to X: ", i)

            self.CNOT(self.qr_X[i], self.qr_left_right[self._get_j(self.A[2], i)])
            self.MeasZ(self.qr_Z[i], cr_z[i])
            #self.Idle(self.qr_right[i])

            self.connectivity_dict_L[i].add(self._get_j(self.A[2], i))

        self.add_barrier()

        #Round 8
        for i in range(int(self.n/2)):
            self.MeasX(self.qr_X[i], cr_x[i])
            self.InitZ(self.qr_Z[i])
            #self.Idle(self.qr_left[i])
            #self.Idle(self.qr_right[i])

        self.add_barrier()

    
    def multiple_syndrome_measurement(self, T):
        """
        Multiple syndrome measurement cycle circuit
        """
        if self.basis == "x":
                self.qc.h(self.qr_left_right)

        for i in range(T):
            self.qc.barrier()
            self.depth_8_syndrome_measurement(i)
        
        self.qc.barrier()
        #final readout
        if self.basis == "x":
                self.qc.h(self.qr_left_right)

        for i in range(self.n):  
            self.qc.measure(self.qr_left_right[i], self.cr_final[i])
            #self.qc.reset(self.qr_left_right[i])


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
        return NotImplementedError
    
    def is_cluster_neutral(self, nodes):
        return NotImplementedError
    
    def string2nodes(self, string, **kwargs):
        return NotImplementedError
    


    def stim_detectors(self):
            """
            Constructs detectors and logicals required for stim.

            Returns:
                detectors (list[dict]): dictionaries containing
                a) 'clbits', the classical bits (register, index) included in the measurement comparisons
                b) 'qubits', the qubits (list of indices) participating in the stabilizer measurements
                c) 'time', measurement round (int) of the earlier measurements in the detector
                d) 'basis', the pauli basis ('x' or 'z') of the stabilizers
                logicals (list[dict]): dictionaries containing
                a) 'clbits', the classical bits (register, index) included in the logical measurement
                b) 'basis', the pauli basis ('x' or 'z') of the logical
            """

            detectors = []
            logicals = []

            ## 0th round of measurements
            if self.basis == "x":
                reg = "round_0_x_bits"
                for stabind, stabilizer in enumerate(self.x_stabilizers):
                    det = {"clbits": []}
                    for gauge_ind in self._gauges4stabilizers[0][stabind]:
                        det["clbits"].append((reg, gauge_ind))
                    det["qubits"] = stabilizer.copy()
                    det["time"] = 0
                    det["basis"] = "x"
                    detectors.append(det)

            else:
                reg = "round_0_z_bits"
                for stabind, stabilizer in enumerate(self.z_stabilizers):
                    det = {"clbits": []}
                    for gauge_ind in self._gauges4stabilizers[1][stabind]:
                        det["clbits"].append((reg, gauge_ind))
                    det["qubits"] = stabilizer.copy()
                    det["time"] = 0
                    det["basis"] = "z"
                    detectors.append(det)

            # adding first x and then z stabilizer comparisons
            for j, basis in enumerate(["x", "z"]):
                for t in range(
                    1, self.T
                ):  # compare stabilizer measurements with previous in each round
                    reg_prev = "round_" + str(t - 1) + "_" + basis + "_bits"
                    reg_t = "round_" + str(t) + "_" + basis + "_bits"
                    for gind, gs in enumerate(self._gauges4stabilizers[j]):
                        det = {"clbits": []}
                        for gauge_ind in gs:
                            det["clbits"].append((reg_t, gauge_ind))
                            det["clbits"].append((reg_prev, gauge_ind))
                        det["qubits"] = self._stabilizers[j][gind].copy()
                        det["time"] = t
                        det["basis"] = basis
                        detectors.append(det)

            ## final measurements
            if self.basis == "x":
                reg_prev = "round_" + str(self.T - 1) + "_x_bits"
                reg_T = "final_readout"
                for stabind, stabilizer in enumerate(self.x_stabilizers):
                    det = {"clbits": []}
                    for q in stabilizer:
                        det["clbits"].append((reg_T, q))
                    for gauge_ind in self._gauges4stabilizers[0][stabind]:
                        det["clbits"].append((reg_prev, gauge_ind))
                    det["qubits"] = stabilizer.copy()
                    det["time"] = self.T
                    det["basis"] = "x"
                    detectors.append(det)
                logicals.append(
                    {
                        "clbits": [(reg_T, q) for q in sorted(self.logical_x[0])],
                        "basis": "x",
                    }
                )
            else:
                reg_prev = "round_" + str(self.T - 1) + "_z_bits"
                reg_T = "final_readout"
                for stabind, stabilizer in enumerate(self.z_stabilizers):
                    det = {"clbits": []}
                    for q in stabilizer:
                        det["clbits"].append((reg_T, q))
                    for gauge_ind in self._gauges4stabilizers[1][stabind]:
                        det["clbits"].append((reg_prev, gauge_ind))
                    det["qubits"] = stabilizer.copy()
                    det["time"] = self.T
                    det["basis"] = "x"
                    detectors.append(det)
                logicals.append(
                    {
                        "clbits": [(reg_T, q) for q in sorted(self.logical_z[0])],
                        "basis": "z",
                    }
                )

            return detectors, logicals



if __name__ == "__main__":
    code = GrossCode()
    circuit = GrossCodeCircuit(code)
    #circuit.qc.draw(output='mpl', filename='gross_code_circuit.png', vertical_compression='high', scale=0.3, fold=500)
    circuit.verify_connectivity()
    pass