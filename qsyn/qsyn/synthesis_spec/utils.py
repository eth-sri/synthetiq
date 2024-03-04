from typing import List, Tuple
import numpy as np
import cirq

if __name__ == "synthesis_spec.utils":
    
    from util.utils import basis_rep
else:
    from qsyn.util.utils import *


class IOpairs():
    def __init__(self, io_pairs: List[Tuple[np.array, np.array]] = None):
        if io_pairs is None:
            self.io_pairs = []
        else:
            self.io_pairs = io_pairs

    def get_io_pairs(self) -> List[Tuple[np.array, np.array]]:
        return self.io_pairs

    def get_identity_input_specs(self):
        collect = list()
        for in_arr, out_arr in self.io_pairs :
            if basis_rep(in_arr) == basis_rep(out_arr):
                collect.append(basis_rep(in_arr)) # save in string

        return collect
    

    def append(self, io_pair: Tuple[np.array, np.array]):
        self.io_pairs.append(io_pair)

    def __str__(self) -> str:
        string_builder = ""
        for idx, io_pair in enumerate(self.io_pairs):
            in_, out_ = io_pair
            string_builder += (
                "(" + str(idx) + ")" + "\n"
                + f"{'-input' :<8}{':':^2}" + cirq.qis.dirac_notation(in_)
                + "\n"
                + f"{'-output':<8}{':':^2}" + cirq.qis.dirac_notation(out_)
                + "\n"
            )
        return string_builder


def is_valid_amplitude(complex_array: np.ndarray):
    norm = np.sum(np.abs(complex_array) ** 2)
    return np.isclose(norm, 1, atol=1e-4)
