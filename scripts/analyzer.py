from qiskit import QuantumCircuit
import os
import numpy as np
from qiskit.tools.visualization import circuit_drawer
from tqdm import tqdm
from qiskit.quantum_info import Operator

OUTPUT_FOLDER = "data/output/"
MONITORING_FOLDER = "monitoring/"
NON_STANDARD_GATES = {
    'scz': ('crz(pi)', Operator(np.array([[1, 0, 0, 0],
                       [0, 1, 0, 0],
                       [0, 0, 1, 0],
                       [0, 0, 0, 1j]]))),
    'U': ('crx(pi)', Operator(np.array([[-0.35355339+0.35355339j,  0.35355339+0.35355339j, 0.35355339+0.35355339j,  0.35355339-0.35355339j],
                                        [ 0.35355339-0.35355339j,  0.35355339+0.35355339j,-0.35355339-0.35355339j,  0.35355339-0.35355339j],
                                        [ 0.35355339-0.35355339j, -0.35355339-0.35355339j, 0.35355339+0.35355339j,  0.35355339-0.35355339j],
                                        [ 0.35355339-0.35355339j,  0.35355339+0.35355339j, 0.35355339+0.35355339j, -0.35355339+0.35355339j]]))),
}

class Circuit:
    def __init__(self, filename) -> None:
        with open(filename, 'r') as file:
            qasm_str = file.read()
        for replace_gate in NON_STANDARD_GATES:
            qasm_str = qasm_str.replace(replace_gate, NON_STANDARD_GATES[replace_gate][0])
        self.circuit = QuantumCircuit.from_qasm_str(qasm_str)
        for index, gate in enumerate(self.circuit.data):
            for replace_gate in NON_STANDARD_GATES:
                if gate[0].name == NON_STANDARD_GATES[replace_gate][0].split('(')[0] and gate[0].params[0] == np.pi:
                    # set gate to the correct gate
                    self.circuit.data[index] = (NON_STANDARD_GATES[replace_gate][1], gate[1], gate[2])

        self.filename = filename
        self.score = float(os.path.basename(filename).split("-")[0])
        self.t_depth = self.circuit.depth(lambda gate: gate[0].name in ['t', 'tdg'])
        self.cx_depth = self.circuit.depth(lambda gate: gate[0].name == 'cx')
        self.cx_count = np.count_nonzero(np.array([el[0].name for el in self.circuit.data]) == "cx")
        gates_names = np.array([el[0].name for el in self.circuit.data])
        self.t_count = np.count_nonzero(gates_names == "t") + np.count_nonzero(gates_names == "tdg")
        self.count = float(os.path.basename(filename).split("-")[1])
        self.gates = len(gates_names)

def main_analysis(circuit_folder):
    t_depth = []
    t_count = []
    gates = []
    best_t_depth_circ = None
    best_t_count_circ = None
    best_cx_depth_circ = None
    for file in tqdm(os.listdir(circuit_folder)):
        circuit = Circuit(os.path.join(circuit_folder,file))
        t_depth.append(circuit.t_depth)
        gates.append(circuit.gates)
        t_count.append(circuit.t_count)

        if best_t_depth_circ is None:
            best_t_depth_circ = circuit
        condition2 = circuit.t_depth < best_t_depth_circ.t_depth
        condition3 = circuit.t_depth == best_t_depth_circ.t_depth and circuit.t_count < best_t_depth_circ.t_count
        condition4 = circuit.t_depth == best_t_depth_circ.t_depth and circuit.t_count == best_t_depth_circ.t_count and circuit.score < best_t_depth_circ.score
        if condition2 or condition3 or condition4:
            best_t_depth_circ = circuit
        
        if best_t_count_circ is None:
            best_t_count_circ = circuit
        condition2 = circuit.t_count < best_t_count_circ.t_count
        condition3 = circuit.t_count == best_t_count_circ.t_count and circuit.t_depth < best_t_count_circ.t_depth
        condition4 = circuit.t_count == best_t_count_circ.t_count and circuit.t_depth == best_t_count_circ.t_depth and circuit.score < best_t_count_circ.score
        if condition2 or condition3 or condition4:
            best_t_count_circ = circuit
            
        # do the same for cx depth
        if best_cx_depth_circ is None:
            best_cx_depth_circ = circuit
        condition2 = circuit.cx_depth < best_cx_depth_circ.cx_depth
        condition3 = circuit.cx_depth == best_cx_depth_circ.cx_depth and circuit.cx_count < best_cx_depth_circ.cx_count
        condition4 = circuit.cx_depth == best_cx_depth_circ.cx_depth and circuit.cx_count == best_cx_depth_circ.cx_count and circuit.score < best_cx_depth_circ.score
        if condition2 or condition3 or condition4:
            best_cx_depth_circ = circuit

    t_depth = np.array(t_depth)
    t_count = np.array(t_count)
    gates = np.array(gates)
    return t_depth, t_count, gates, best_t_count_circ, best_t_depth_circ, best_cx_depth_circ

def analyse_circuits(circuit_folder, save_folder=None):
    if save_folder is None:
        save_folder = os.path.join(MONITORING_FOLDER, OUTPUT_FOLDER, circuit_folder)
    os.makedirs(save_folder, exist_ok=True)

    t_depth, t_count, gates, best_t_count_circ, best_t_depth_circ, best_cx_depth_circ = main_analysis(OUTPUT_FOLDER + circuit_folder)

    with open(os.path.join(save_folder, "analysis.txt"), "w") as f:
        print("---- Analysis of circuits found ----", file=f)

        print(f"Min gate count: {np.min(gates)}", file=f)
        print(f"CX Depth: {best_cx_depth_circ.cx_depth}", file=f)
        
        print("-" * 40, file=f)
        print(f"Best found T-depth: {min(t_depth)}", file=f)
        print(f"Filename: {best_t_depth_circ.filename}", file=f)

        print("-" * 40, file=f)
        print(f"Best found T-count: {min(t_count)}", file=f)
        print(f"Filename: {best_t_count_circ.filename}", file=f)
        print("-" * 40, file=f)

    best_t_count_circ.circuit.qasm(formatted=True, filename=os.path.join(save_folder, "tcount.qasm"))
    best_t_depth_circ.circuit.qasm(formatted=True, filename=os.path.join(save_folder, "tdepth.qasm"))
    best_cx_depth_circ.circuit.qasm(formatted=True, filename=os.path.join(save_folder, "cxdepth.qasm"))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        prog="Circuit Analysis",
        description="Analyzes the output of Synthetiq in a given folder"
    )
    parser.add_argument("circuit_folder")
    parser.add_argument("--save-folder", type=str, help="Where to store the best circuit", default=None)
    args = parser.parse_args()
    analyse_circuits(
        circuit_folder=args.circuit_folder,
        save_folder=args.save_folder,
    )
