import numpy as np


class Matrix(object):

    def __init__(self, values, coeff_modulus=0):
        self.values = np.array(values)
        self.n_rows = self.values.shape[0]
        self.n_cols = self.values.shape[1]
        self.coeff_modulus = coeff_modulus
        if coeff_modulus != 0:
            self.values %= self.coeff_modulus

    def __repr__(self):
        return str(self.values)

    def kernel_image(self):
        B, Q, Qi, k = self.row_echelon()
        return Qi.T[:, k:], B.T[:, :k]

    def row_echelon(self):
        B = self.values.copy()
        Q = np.identity(self.n_rows)
        Qi = np.identity(self.n_rows)
        k = 0
        l = 0
        while k < self.n_rows:
            while l < self.n_cols and not (np.any(B[k:, l])):
                l += 1
            if l == self.n_cols:
                break
            B, Q, Qi = self._row_reduce(B, Q, Qi, k, l)
            k += 1
        return B, Q, Qi, k

    def is_square(self):
        return self.n_rows == self.n_cols

    def _row_reduce(self, B, Q, Qi, k, l):
        while np.any(B[k + 1:, l]):
            B, Q, Qi = self._row_prepare(B, Q, Qi, k, l)
            B, Q, Qi = self._partial_row_reduce(B, Q, Qi, k, l)
        return B, Q, Qi

    def _row_prepare(self, B, Q, Qi, k, l):
        (a, i) = self._smallest_nonzero_index(B[:, l], k)
        B[[i, k], :] = B[[k, i], :]
        Qi[[i, k], :] = Qi[[k, i], :]  # row swap
        Q[:, [i, k]] = Q[:, [k, i]]  # column swap
        return B, Q, Qi

    @staticmethod
    def _smallest_nonzero_index(v, k):
        # TODO: This pattern is bad
        try:
            alpha = min(abs(v[k:][np.nonzero(v[k:])]))
            i = min(i for i in range(k, len(v)) if abs(v[i]) == alpha)
            return alpha, i
        except:
            return 0, np.nan  # minNonZero causes this sometimes

    def _partial_row_reduce(self, B, Q, Qi, k, l):
        for i in range(k + 1, Q.shape[0]):
            q = (B[i, l] // B[k, l])
            if self.coeff_modulus != 0:
                q %= self.coeff_modulus
            # row add i,k,-q
            B[i] += (-q * B[k])
            Qi[i] += (-q * Qi[k])  # row add
            Q[:, k] += (q * Q[:, i])  # column add (note i,k are switched)
            if self.coeff_modulus != 0:
                B[i] %= self.coeff_modulus
                Qi[i] %= self.coeff_modulus
                Q[:, k] %= self.coeff_modulus
        return B, Q, Qi
