# Introduction
This repository contains the code for the paper 'Synthetiq: Fast and Versatile Quantum Circuit Synthesis'. Synthetiq is a tool to synthesize quantum circuits implementing a given (partial) specification over arbitrary finite gate sets and is faster and more versatile than existing works. This README is a guide on installing and using Synthetiq, as well as replicating the paper's results.

For easy reference, all quantum circuits found by Synthetiq and mentioned in the paper can be found in [`data/circuits`](data/circuits), under the subfolder corresponding to the subsection they are mentioned in. All circuits are in OpenQASM 2.0 format. 

We describe how to reproduce all claims made in the paper. We recall them here from Section 6 in the paper, and indicate for each the relevant figures.
1. The circuit cost rewrite, permutation of qubits in cost and the simplification pass increase Synthetiq performance (Table 2) and our hyperparameters values are optimal on our benchmark (Table 3).
2. Synthetiq finds better implementations than the currently best known ones for RCCCX, CCT, CCiSWAP, C√SWAP, and C√iSWAP (Table 4).
3. Synthetiq can efficiently synthesize circuits with user-supplied custom gate sets, and outperforms the state-of-the-art [Kang and Oh 2023] in 50% of the cases (including 27% of cases where Kang and Oh [2023] fails to return any decomposition) while being equally optimal in all other cases (Figure 3 and Table 5).
4. For completely specified operators over the Clifford+T gate set, Synthetiq outperforms the specialized state-of-the-art [Gheorghiu et al. 2022a; Mosca and Mukhopadhyay 2021]. Synthetiq is able to find circuits for more operators, those circuits are often more efficient and typically found one to two orders of magnitude faster (Table 6 and Figure 4).
5. For approximate synthesis over the Clifford+ T gate set, Synthetiq is 6 times faster than the state-of-the-art approach specialized to this task [Gheorghiu et al. 2022b] for complex multi-qubit operators and, while slower, can find circuits that are on par with Gheorghiu et al. [2022b] for simpler single qubit operators (Table 7).
6. Using Synthetiq to synthesize small components of a bigger circuit allows for more efficient implementations. Specifically, by synthesizing a relative phase carry operator, we can reduce the T-count of the Cirq Adder [Cirq 2023] by more than 3x (Figure 6).

For all claims, we provide all the benchmarks we used. Further, for Claim 3, we provide the code for our baseline [Kang and Oh 2023] as well as instructions on how to run it. Finally, for Claims 4 and 5, we cannot provide the code for the baselines, as it was communicated to us privately and the authors declined to make it public. We however provide all their raw results on our benchmarks.

# Hardware Dependencies
No specific hardware is required to run Synthetiq. However, re-running our complete evaluation requires more than a week of compute on a machine with 64 cores.
To allow evaluation of this artifact, we therefore provide (i) instructions to re-run it completely, (ii) instructions to re-run only a small subset of it (with only a few hours of compute) and (iii) the raw results of a complete re-run, as well as instructions on how to reproduce our paper's figures and tables from those. Note that as those raw results are from a re-run and Synthetiq is based on a random process, we expect to see slight variations in the final results.

# Getting Started Guide 
## Installation
Synthetiq can either be installed directly or through Docker. After unzipping the artifact, run the following from its main folder.

### Direct Installation
Synthetiq requires gcc, make, and Python for direct installation. The following has been tested on Ubuntu 20 and 22. 

Install gcc and make with:
```bash
sudo apt-get update
sudo apt-get install -y gcc-11 make
```
Python installation via [Conda](https://docs.conda.io/projects/miniconda/en/latest/) is recommended. We note that Python is only necessary to run post-processing. It is not necessary to run Synthetiq itself. 

Follow the Conda installation instructions and execute:
 ```bash
 conda create -n synthetiq python=3.10
 conda activate synthetiq
 python -m pip install -r requirements.txt
 ```
Then, install Synthetiq:
```bash
make clean
make all
```

### Installation Using Docker
First, download and install [Docker](https://www.docker.com/get-started/) per the instructions. Then, build and run Synthetiq in Docker:
```bash
docker build -t synthetiq .
docker run -it -p 8888:8888 synthetiq /bin/bash
```

## Testing the Installation
Synthetiq uses a command-line interface with one mandatory argument, the operator specification it should build an implementation for.

For instance, running 
```bash
  ./bin/main cx.txt
```
will run Synthetiq for the specification in the file [`data/input/cx.txt`](data/input/cx.txt), using the Clifford+T gate set. This run will last either 100s or until 10 circuits implementing this specification are found. You should get an output similar to the following:

```bash
Input file: data/input/cx.txt
Output folder: data/output/cx/
Format: txt
Gate set: data/gates/CliffordT
Composite gate folder: data/gates/composite_gates
Number of threads: 1
Number of circuits to find: 10
Number of ancillas: 0
Number of dirty qubits: 0
Optimal t-count: -1
Optimal t-depth: -1
Optimal gate count: -1
Cost required: -1
Time allowed: 100
Start temp base: 0.1
Epsilon: 1e-06
Beta: 0.05
Min factor: 2.5
Max factor: 3.5
Pcomp: 0.2
Pid: 0.3
Enable permutations: 1
Do resynth: 1
Simple cost: 0
N norm: 80
Iterations factor: 40
Update gate scheme: 1
Save times: 0
Expand composite: 0
Qubit independent: 1
Save all circuits: 0
Save any circuit: 1
Times file: data/times.csv
Optimization number: 12
Optimize depth: 1
Inverse independent: 0
Depth gates: t tdg 
data/output/cx/ 0.0303801 1
Gates: 57 -> 23
T-count: 16 -> 4
T-depth: 12 -> 3
data/output/cx/ 0.0134407 2
Gates: 36 -> 1
T-count: 12 -> 0
T-depth: 10 -> 0
data/output/cx/ 6.3587e-05 3
Gates: 1 -> 1
T-count: 0 -> 0
T-depth: 0 -> 0
data/output/cx/ 4.5299e-05 4
Gates: 3 -> 1
T-count: 0 -> 0
T-depth: 0 -> 0
data/output/cx/ 3.732e-05 5
Gates: 3 -> 1
T-count: 0 -> 0
T-depth: 0 -> 0
data/output/cx/ 4.0085e-05 6
Gates: 3 -> 1
T-count: 0 -> 0
T-depth: 0 -> 0
data/output/cx/ 6.5157e-05 7
Gates: 1 -> 1
T-count: 0 -> 0
T-depth: 0 -> 0
data/output/cx/ 8.5035e-05 8
Gates: 3 -> 1
T-count: 2 -> 0
T-depth: 2 -> 0
data/output/cx/ 4.1324e-05 9
Gates: 1 -> 1
T-count: 0 -> 0
T-depth: 0 -> 0
data/output/cx/ 5.571e-05 10
Gates: 3 -> 1
T-count: 0 -> 0
T-depth: 0 -> 0
```

Further, all found implementations will be stored in the folder [`data/output/cx`](data/output/cx). Each implementation is returned in the OpenQASM 2.0 format, and has the following naming convention: `[Total Cost]-[T-count]-[T-depth]-[Id-1]-[Id-2].qasm`. There should be 10 such files.

To check your python installation, you can run
```bash
  python scripts/main.py
```
this will run the complete post-processing of the results. The output should show progress bars along with some circuits printed in OpenQASM format.

# Step-by-Step Instructions to Reproduce Evaluation
In the following, we describe how to reproduce the evaluation results from our paper, "Synthetiq: Fast and Versatile Quantum Circuit Synthesis". 

As re-running all of our experiments takes more than a week, we provide multiple options for reproducibility:
1. Using our provided raw results and only run post-processing, which we expect to take a few seconds. To do so, skip to section [Processing of Raw Results](#processing-of-raw-results). We note that the provided raw results are from a rerun on our machine. Therefore, the slight variations in time are the ones expected based on the random variation between runs.
2. Running only a subset of our experiments, which we found to take 2h30min on a laptop (when setting threads=4). To do so, follow section [Re-Running Small Experiments](#re-running-small-experiments).
3. Running all experiments. To do so, follow section [Re-Running All Experiments](#re-running-all-experiments).

## Re-Running All Experiments
Before running anything in this section, we recommend deleting all the raw results stored in [`data/output`](data/output), [`data/times.csv`](data/times.csv) and [`data/optimization.csv`](data/optimization.csv).
We provide a single script to run the complete evaluation, producing all raw results used in the paper. It can be run using 
```bash
bash scripts/main.sh
```
in the folder where Synthetiq is installed. This script assumes you have 64 cores available. 
If this is not the case, make sure to update this number at the top of the file. Note that this may yield different results than those presented in the paper. If you only want to run a specific experiment, you can simply comment out the lines of the script corresponding to other experiments. 

Further, for Claim 3, it is necessary to also run our baseline [Kang and Oh 2023]. To do so, delete [`qsyn/output`](qsyn/output) and run the following commands from the [`qsyn`](qsyn) subfolder:
```bash
docker build -t qsyn .
bash main.sh
```
Note that this `main.sh` assumes you have 64 cores available. If this is not the case, either retain the raw output of the baseline we provide or appropriately adjust the file to use less cores. The gate counts of this baseline cannot be extracted automatically, since the found circuits are only printed in a difficult-to-parse format. After the runs are finished, you have to manually inspect the files in the `output` folder and count the number of gates at the end of the raw output files. Put those counts in the `kang_counts` variable at the top in the notebook [`notebooks/post_processing/63.ipynb`](notebooks/post_processing/63.ipynb).

## Re-Running Small Experiments
As above, before running anything in this section, we recommend deleting all the raw results stored in [`data/output`](data/output), [`data/times.csv`](data/times.csv) and [`data/optimization.csv`](data/optimization.csv).

We provide a single script to run this shorter evaluation, producing raw results for Sections 6.3 (only Table 5), 6.4 (only the CCX line in Table 4), 6.5 (only the results which took Synthetiq less than 1s to find) and 6.6 in the paper. It can be run using 
```bash
bash scripts/main_small.sh
```
in the folder where Synthetiq is installed. This script assumes you have 64 cores available. 
If this is not the case, make sure to update this number at the top of the file. Note that this may yield different results than those presented in the paper. If you only want to run a specific experiment, you can simply comment out the lines of the script corresponding to other experiments. 

After running this script, you can run the following steps of the processing (other steps require more data, that is only produced by running all experiements):

### Section 6.3: Mode: Custom Gates
In the notebook [`notebooks/post_processing/63.ipynb`](notebooks/post_processing/63.ipynb), you should run the first 4 cells and the last one. This will produce Table 5 in [`data/paper/63/table5.csv`](data/paper/63/table5.csv) and save the best circuits in [`data/paper/63/stackexchange`](data/paper/63/stackexchange/).

### Section 6.4: Mode: Clifford+T Gate Set
In the notebook [`notebooks/post_processing/64.ipynb`](notebooks/post_processing/64.ipynb) you can run the first five cells, which will produce Table 6 in [`data/paper/64/comparison/table6.csv`](data/paper/64/comparison/table6.csv) and save the best circuit found for CCX in [`data/paper/64/comparison`](data/paper/64/comparison/).

### Section 6.5: Mode: Approximate Circuit Synthesis
Running the whole post-processing in [`notebooks/post_processing/65.ipynb`](notebooks/post_processing/65.ipynb) produces Table 7 in [`data/paper/65/table7.csv`](data/paper/65/table7.csv) and save the best circuits in [`data/paper/65`](data/paper/65/).


### Section 6.6: Mode: Relative Phase Operators

Running the whole post-processing in [`notebooks/post_processing/66.ipynb`](notebooks/post_processing/66.ipynb) produces Figure 6 in [`data/paper/fig6`](data/paper/fig6/).


## Processing of Raw Results
It is possible to run all the post-processing at once by running 
```bash
python scripts/main.py
```
This will reproduce all the results from our paper, and put them in [`data/paper`](data/paper), ordered by section. You can there read and check the reproduced results. 

You can also run this post-processing step by step, using the jupyter notebooks we provide.

If you use the direct installation, you can run the notebooks using:
```bash
jupyter notebook
```

If running through docker, you can run the notebooks using:
```bash
jupyter notebook --ip 0.0.0.0 --no-browser --allow-root
```
, clicking one of the links shown and navigating to the relevant notebook.

### Section 6.1: Hyperparameter Optimization
All the necessary post-processing to reproduce Table 2 and results as shown in Table 3 is shown in [`notebooks/post_processing/61.ipynb`](notebooks/post_processing/61.ipynb).

### Section 6.2: Better Operator Decompositions
All the necessary post-processing to produce Table 4 is shown in [`notebooks/post_processing/62.ipynb`](notebooks/post_processing/62.ipynb). As mentioned in the introduction, we also provide in [`data/circuits/62`](data/circuits/62) the best decomposition found for each circuit by Synthetiq as a `.qasm` file.

### Section 6.3: Mode: Custom Gates
You can follow the notebook [`notebooks/post_processing/63.ipynb`](notebooks/post_processing/63.ipynb) to produce Figure 3 and Table 5.

As above, we also provide in [`data/circuits/63`](data/circuits/63) the best decomposition found by Synthetiq for each of the stackexchange problems.

### Section 6.4: Mode: Clifford+T Gate Set
As mentioned above, the code of the baselines we use for this experiment (Mosca and Mukhopadhyay [2021] and Gheorghiu et al. [2022a]) was sent to us in private communication, and we did not get permission from the authors to include their code here. We therefore provide the raw outputs of these runs in [`data/baselines`](data/baselines).

The post-processing in [`notebooks/post_processing/64.ipynb`](notebooks/post_processing/64.ipynb) then produces Figure 4 and Table 6.

Again, we also provide in [`data/circuits/64`](data/circuits/64) the best decomposition found by Synthetiq for each of the operators in Table 6.

### Section 6.5: Mode: Approximate Circuit Synthesis
As above, we could not include the code of the baseline Gheorghiu et al. [2022b], and therefore provide its raw output in [`data/baselines`](data/baselines).

The post-processing in [`notebooks/post_processing/65.ipynb`](notebooks/post_processing/65.ipynb) then produces Table 7.

Again, we also provide in [`data/circuits/65`](data/circuits/65) the best decomposition found by Synthetiq for each of the operators in Table 7.

### Section 6.6: Mode: Relative Phase Operators

The post-processing in [`notebooks/post_processing/66.ipynb`](notebooks/post_processing/66.ipynb) produces Figure 7.

### Notes
1. The output of Synthetiq is based on a random process. We use seeds to help with reproducibility, but requiring the different threads to synchronize enough for appropriate seeding proved terrible for performance. Therefore, our results are only properly seeded when using only one thread. If you use more than one thread, variation between runs can lead to minor differences between them. 
2. For the hyperparameter optimization phase, we rely on several randomly generated operators. You can regenerate these operators by deleting the [`data/input/61`](data/input/61) files and running
```bash
./bin/comparison_generator
```
3. You can regenerate the 30 3-qubit permutations we evaluate on in Section 6.2 of the paper using [`notebooks/clifford_equivalence.ipynb`](notebooks/clifford_equivalence.ipynb). In the code, you can find how we ensured that all other permutations can be transformed into one of these 30 permutations using Clifford operations.

# Reusability Guide 

## Using Synthetiq
Synthetiq's command-line interface has one mandatory argument, the operator specification it should build an implementation for.

For instance, running 
```bash
  ./bin/main cx.txt
```
will run Synthetiq for the specification in the file [`data/input/cx.txt`](data/input/cx.txt), using the Clifford+T gate set. This run will last either 100s or until 10 circuits implementing this specification are found. Once a circuit is found, you will see some basic information about it in the output in the following format:
```
Gates: [Number of gates before simplification] -> [Number of gates after simplification]
T-count: [T-count before simplification] -> [T-count after simplification]
T-depth: [T-depth before simplification] -> [T-depth after simplification]
```
All found implementations are stored in the folder [`data/output/cx`](data/output/cx). Each implementation is returned in the OpenQASM 2.0 format, and has the following naming convention: `[Total Cost]-[T-count]-[T-depth]-[Id-1]-[Id-2].qasm`. Here, Id-1 and Id-2 are numbers that can be safely ignored when interpreting the output. If no circuit satisfying the specification are found, no such files are created.

We now list the main optional arguments Synthetiq provides, along with their default value:

| Option | Default Value | Explanation |
|-----------------|-----------------|-----------------|
| --output / -o | Mirrored from input | Specifies the output folder in which to place the found circuits, relative to the [`data/output`](data/output) folder |
| --threads / -h | 1 | Number of threads to use for parallellization |
| --time / -t | 100 | Time in seconds Synthetiq is allowed to run for |
| --circuits / -c | 10 | Number of implementations to find before exiting |
| --gate-set / -gs | CliffordT | The gate set folder to use, relative to the [`data/gates`](data/gates) folder |
| --composite-gates / -cg | composite_gates | The composite gate folder to use, relative to the [`data/gates`](data/gates) folder. The `composite_gates` folder is empty by default. | 
| --depth-gates | t,tdg | List of gates for which to optimize the depth in the simplification pass |
| --epsilon / -eps | 1e-6 | Float indicating to what precision the implementation should match the specification |
| --gates / -g | None | If provided, Synthetiq will assume a circuit with the specified number of gates has been found (this is used for setting the circuit size bounds) |
| --ancilla / -a | 0 | Number of ancillae Synthetiq has to use |
| --dirty / -d | 0 | Number of dirty qubits Synthetiq has to use |
| --expand / -e | boolean (and by default not set) | When set, composite gates will be expanded before running the simplification pass |
| --tcount / -tc | None | Synthetiq will discard and ignore any implementation with a higher T-count |
| --tdepth / -td | None | Synthetiq will discard and ignore any implementation with a higher T-depth |
| --gate-count / -gc | None | Synthetiq will discard and ignore any implementation with a higher gate count |
| --cost-required / -cr | None | Synthetiq will discard and ignore any implementation with a higher total cost |
| --absolute-input | false | The input file specified is an absolute path instead of a path relative to data/input |
| --absolute-output | false | The output folder specified is an absolute path instead of a path relative to data/output |
| --absolute-gates | false | The gate folder specified is an absolute path instead of a path relative to data/gates |


For instance, setting more arguments explicitly for the example above results in the following command:
```bash
./bin/main cx.txt --output cx --threads 1 --time 100 --circuits 10 --gate-set CliffordT -cg composite_gates
```

### Writing operators specifications
The first argument of Synthetiq is the operator specification, as a txt file. For this file, we use the following human-readable format:
```
[matrix name]
[number of qubits]
[Operator specification (U in our paper), complex matrix]
[Boolean Cover (M in our paper), boolean matrix]
```
For instance, the following operator acts on two qubits and is only partially specified:
```
example
2
(0,0) (0.5,0) (0,0) (0,0)
(0.5,0) (-0.5,0) (0.5,0) (-0.5,0)
(0.5,0) (-0.5,0) (-0.5,0) (0.5,0)
(0.5,0) (0.5,0) (0,0) (-0.5,0)
0 1 0 0
1 1 1 1
1 1 1 1
1 1 0 1
```
As mentioned above, all input specifications should be placed in the [`data/input`](data/input) folder. You can add subfolders within this folder to structure your workflow. 

**Note**: You can also specify your input operator as a circuit in the OpenQASM-2.0 format. Make sure to use the appropriate extension name `.qasm`. Synthetiq will then automatically generate the matrix of the given circuit and try to find a more efficient implementation than the one given. 

In the [`data/input`](data/input) folder, we provide the specifications of all the operators we used in our paper. 

### Specifying a gate set

For the gate set, we use a very similar syntax. Each gate is specified as a `.txt` file structured as follows:
```
[gate name]
[number of qubits]
[Gate Cost]
[Space separated list of qubits on which the gate acts]
[Gate specification, complex matrix]
```
For the CX gate, this becomes:
```
cx
2
0.1
1 0
(1,0) (0,0) (0,0) (0,0)
(0,0) (1,0) (0,0) (0,0)
(0,0) (0,0) (0,0) (1,0)
(0,0) (0,0) (1,0) (0,0)
```
Your gate set should then be a subfolder in the [`data/gates`](data/gates) folder, and each file within your subfolder should specify a gate in the gate set. All gate sets used in our paper are already available there.

### Specifying composite gates

To use composite gates in the decomposition, you can create another subfolder in the [`data/gates`](data/gates) folder. You can then describe each composite gate as a `.qasm` file, using the OpenQASM 2.0 format to specify its decomposition in the gates from the gate set you use. 

If you are using Docker, make sure to rerun 
```bash
docker build -t synthetiq .
```
if you have added your file to the specified folders. This will update the folders copied to the docker image.

### Running the simplification pass
Synthetiq implements a simplification pass, which is run on all circuits found. It is also possible to run this simplification pass on any circuit in OpenQASM format. For instance, the following command:

```bash
./bin/main_resynth example_circuit_to_simplify.qasm
```
will run the simplification pass on the circuit specified in the file [`data/input/example_circuit_to_simplify.qasm`](data/input/example_circuit_to_simplify.qasm). The output of the command will indicate whether the new T-count and T-depth is better than the one before the simplification. Specifically, the output will look like
```
Gates: [Number of gates before simplification] -> [Number of gates after simplification]
T-count: [T-count before simplification] -> [T-count after simplification]
T-depth: [T-depth before simplification] -> [T-depth after simplification]
```
## Code Structure
All the implementation of Synthetiq is in the folder [`synthetiq`](synthetiq/). We briefly present its main files below:
1. [algo.cpp](synthetiq/algo.cpp) is the entrance point. It parses all arguments and creates all the necessary objects that will be used in a run of Synthetiq, and runs iterations of [mcmc_sa.cpp](synthetiq/mcmc_sa.cpp).
1. [mcmc_sa.cpp](synthetiq/mcmc_sa.cpp) implements the main algorithm of Synthetiq, Simulated Annealing.
1. [circuit.cpp](synthetiq/circuit.cpp) implements our representation of circuits. This implementation relies on gates, implemented in [gate.cpp](synthetiq/gate.cpp). We use [matrix_computer.cpp](synthetiq/matrix_computer.cpp) to efficiently compute this circuits operator as an `Eigen::MatrixXcd` matrix.
1. [partialMatrix.cpp](synthetiq/partialMatrix.cpp) implements our representation of partial specifications.
1. [resynthesis.cpp](synthetiq/resynthesis.cpp) implements our resynthesis algorithm, which is run as part of [algo.cpp](synthetiq/algo.cpp), but can also be run directly from [main_resynth.cpp](synthetiq/main_resynth.cpp).
1. [cost.cpp](synthetiq/cost.cpp) implements the cost functions described in our paper.

## Cite

```
@inproceedings{paradis2024synthetiq,
  title={Synthetiq: Fast and Versatile Quantum Circuit Synthesis},
  author={Paradis, Anouk and Dekoninck, Jasper and Bichsel, Benjamin and Vechev, Martin},
  journal      = {Proc. {ACM} Program. Lang.},
  volume       = {8},
  number       = {{OOPSLA1}},
  year={2024},
}
```