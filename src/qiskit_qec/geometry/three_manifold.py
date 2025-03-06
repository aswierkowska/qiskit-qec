from qiskit_qec.geometry.manifold import Manifold


class ThreeManifold(Manifold):
    """Like Manifold but better"""

    def __init__(self):
        """Init Manifold"""
        super().__init__(dim=3)