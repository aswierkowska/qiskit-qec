from typing import Optional
from qiskit_qec.codes.codebuilders.builder import Builder

class GrossCodeBuilder(Builder):

    def __init__(self, 
                code_length: Optional[int] = 144,
                k: Optional[int] = 12,
                distance: Optional[int] = 12
    ) -> None:
    
        """Initializes a Gross code builder
    
        Args:
            code_length: length of the code
            k: number of logical qubits
            distance: distance of the code
            
        """
        

