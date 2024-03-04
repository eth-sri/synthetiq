import cirq
from typing import List, Union, Optional
from .utils import *
import numpy as np


class InvalidIOException(Exception):
    def __init__(self, io_pairs):
        msg_builder = ""
        invalid_arr_type = None
        invalid_arr = None
        invalid_arr_idx= None
        for pair_idx, io_pair in enumerate(io_pairs):
                    in_nparr, out_nparr = io_pair
                    if not is_valid_amplitude(in_nparr):
                        invalid_arr_type = "input"
                        invalid_arr_idx = pair_idx
                        invalid_arr =  in_nparr
                        break
                    if not is_valid_amplitude(out_nparr):
                        invalid_arr_type = "output"
                        invalid_arr_idx = pair_idx
                        invalid_arr =  out_nparr
                        break 
        super().__init__(f'\n\n{invalid_arr_idx}-th {invalid_arr_type} state vector {invalid_arr} is invalid' 
                         +'\nThe state vector should be normalized.')

    

class Spec():
    FULL_SYN_SPEC = "Full Spec"
    PART_SYN_SPEC = "Partial Spec"

    def __init__(self, ID: str,  target_object: Union[cirq.Circuit, np.array, IOpairs], qreg_size: int, max_inst_number: int, components: List[cirq.Gate], equiv_phase: bool):
        self.ID = ID
        self._target_object = target_object
        # target object must be full-unitary matrix of qreg_size \tiems qreg_size or cirq_object?
        # or IO pair of (|x> \to |x'>) where each |x>,|x'> must be of 2^(qreg_size) vector
        # the pair of (I,O) maybe presented in np_array w.r.t computational basis.
        self._qreg_size = qreg_size
        self.max_inst_number = max_inst_number
        self._components = components
        self._equiv_phase = equiv_phase
        assert isinstance(self._equiv_phase, bool)
        if self.type_of_spec == Spec.PART_SYN_SPEC:
            self.num_of_ios = len(self._target_object.get_io_pairs())
        if self.validate_spec() == False:
            io_pairs = self.spec_object.get_io_pairs().copy()
            raise InvalidIOException(io_pairs)

    @property
    def id(self):
        return self.ID

    def validate_spec(self) -> bool:
        if self.type_of_spec == self.FULL_SYN_SPEC:
            if isinstance(self.spec_object, np.ndarray):
                return cirq.linalg.is_unitary(self.spec_object)
            elif isinstance(self.spec_object, cirq.Circuit):
                return True
        elif self.type_of_spec == self.PART_SYN_SPEC:
            if isinstance(self.spec_object, IOpairs):
                for io_pair in self.spec_object.get_io_pairs():
                    in_nparr, out_nparr = io_pair
                    if not is_valid_amplitude(in_nparr):
                        return False
                    if not is_valid_amplitude(out_nparr):
                        return False
                return True

    @property
    def spec_object(self):
        return self._target_object

    @property
    def qreg_size(self) -> int:
        return self._qreg_size

    def get_max_inst_number(self) -> int:
        return self.max_inst_number

    @ property
    def component_gates(self):
        return self._components

    @ property
    def equiv_phase(self):
        return self._equiv_phase

    @ property
    def type_of_spec(self) -> str:
        if isinstance(self.spec_object, cirq.Circuit) or isinstance(self.spec_object, np.ndarray):
            return self.FULL_SYN_SPEC
        elif isinstance(self.spec_object, IOpairs):
            return self.PART_SYN_SPEC

        raise Exception(
            'Invalid Spec Type due to invalid instance of spec object')

    def print_components(self) -> str:
        return ",".join([str(gate) for gate in self.component_gates])

    def print_spec_object(self) -> str:
        if self.type_of_spec == self.PART_SYN_SPEC:
            return str(self.spec_object)
        elif self.type_of_spec == self.FULL_SYN_SPEC:
            return str(self.spec_object)

    def __str__(self) -> str:
        string_builder = (
            f"{'ID':<16}{self.ID:>20}" + "\n"
            + f"{'QREG SIZE':<16}{self.qreg_size:>20}" + "\n"
            + f"{'MAX INST ALLOWED':<16}{self.max_inst_number:>20}" + "\n"
            + f"{'PHASE EQUIV':<16}" + str(self.equiv_phase) + "\n"
            + f"{'COMPONENTS':<16}" + "\n" + self.print_components() + "\n"
            + f"{'SPEC OBJECT':<16}" + "\n" + self.print_spec_object() + "\n"
        )

        return string_builder
    def __len__(self) -> str :
        if isinstance(self.spec_object, IOpairs):
            return len(self.spec_object.get_io_pairs())
        else : return 0 
