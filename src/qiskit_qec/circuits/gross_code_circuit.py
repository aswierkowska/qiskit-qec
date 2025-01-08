"""Generates circuits for the gross code."""
from qiskit_qec.circuits.code_circuit import CodeCircuit

class GrossCodeCircuit(CodeCircuit):
    
    #TODO: Figure out what inputs are needed
    def __init__(self):
        super().__init__()
    pass

    #prpoably we need these functions
    def x(self):
        """Apply a logical X gate to the code."""
        pass

    def z(self):
        """Apply a logical Z gate to the code."""
        pass

    def syndrome_measurement(self):
        """
        One sydrome measurement round
        """
        pass

    def readout(self):
        """
        Readout of all logical qubits
        """
        pass