import numpy as np

class BBus:
    def __init__(self, DBAR, DLIN, ref_bar):
        self.DLIN = DLIN
        self.b_ = np.zeros((len(DBAR), len(DBAR)))

        self.create_triangular_admittance()
        self.create_main_diagonal()

        self.b_[ref_bar, ref_bar] = np.inf

    def create_triangular_admittance(self):
        for i, row in enumerate(self.DLIN):
            _b_k = int(row['DE'] - 1)
            _b_m = int(row['PARA'] - 1)
            self.b_[_b_k, _b_m] += row['SUSCEPTANCIA']
            self.b_[_b_m, _b_k] += row['SUSCEPTANCIA']

    def create_main_diagonal(self):
        for k in range(self.b_.shape[0]):
            self.b_[k, k] = sum(self.b_[k, :]) * (-1)