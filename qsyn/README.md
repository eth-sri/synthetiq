

# Modular Component-based Quantum Circuit Synthesis

This is repository for implementation of paper "Modular Component-based Quantum Circuit Synthesis".
Our algorithm produce synthesis of quantum circuit, given by user-provided in/output spec and component gates.

**For published-version of artifact and artifact manual**, see zenodo repository : https://doi.org/10.5281/zenodo.7710436


## Setup

We present two option for setup : (i) by building Docker image or (ii) by installing required python libraries via `pip`. Among two options, we recommend using Docker so that user's local environment is separated. 


### Via Docker
We packaged our artifact in Docker.
Following command will build the Dockerfile, in few minutes (<10 min).

```
docker build -t qsyn -f Dockerfile .
```


Note that `-t` option denotes tag name for build image.
After the build, run docker image with following command:

```
docker run -ti qsyn
```

### Via `pip` install


Environment for our program also can be prepared, only by installing required python libraries listed in `requirements.txt` via `pip`, instead using Docker.
Installation will be done by following command :

```
pip3 install -r requirements.txt
```



## Verifying Install

If the environment is properly prepared, under entry directory(`/workdir/` if using docker) we will be able to run following command (which synthesizes our working example `GHZ_from_100` appeared Section 3 in our paper).

```
python3 qsyn/run_single.py --benchmark GHZ_from_100  --mode Ours
```




### Running Single Synthesis Problem


Each run of single benchmark is done by following command :

```
python3 qsyn/run_single.py --benchmark [Benchmark]  --mode [Mode]
```

Two option should be specified in the command:
- For `[Benchmark]`, insert name of one of `.json` files under `/benchmarks/` folder, which is one of following items
```
three_superpose , M_valued , GHZ_from_100 , GHZ_by_iSWAP, GHZ_by_QFT , GHZ_Game , W_orthog , W_phased , W_four ,  cluster , bit_measure, flip 	, teleportation 	, indexed_bell 	,toffoli_by_sqrt_X	, QFT , draper
```
Note that name is aligned with benchmarks in Table~3.
	
- For `[Mode]`, insert name of  synthesis algorithm mode. The mode is one of following four (refer Section~5.1 of our paper for more detail) :
	- `Ours` :  Main setting of our modular synthesis
	- `Ours_no_prune` : Ours with no module-level pruning.
	- `Base`  : Naive-BFS-based synthesis algorithm with basic pruning.
	-  `Base_no_prune` : `Base`  wihtout pruning. This setting is included to check that  `Base` is not very weak.

For example, following command runs synthesis problem `GHZ_from_100` with our main algorithm mode `Ours`.

```
python3 qsyn/run_single.py --benchmark GHZ_from_100  --mode Ours 
```


### Reproducing Table~3 (in our paepr)

Whole benchmark run for reproduction of Table~3 is done by running `run_all.py` file.
Following command will run whole benchmarks under `/benchmarks` for each mode of synthesis:
```
python3 qsyn/run_all.py  --mode [Mode] 
```
where `[Mode]` specifies setting for synthesis, as one of `Ours`, `Ours_no_prune`, `Base`, `Base_no_prune`.

### Giving New Synthesis Problem

In this section, we explain how to prepare (new) synthesis problem that our program will adopt.
Users need to specify several information of synthesis problem, including input/output state vectors and component gates in `.json` format.
Following is template for benchmark `.json` file.
```
{
	"ID": [Name of Synthesis Problem],
	"Spec Type": "Partial Spec",
	"Spec": {
		"qreg_size": [Size of qregister],
		"max_inst_allowed": 10,
		"spec_object": [
			{
			"input":  [input state vector],
			"output": [output state vector]
			},
			{
			"input":  [input state vector],
			"output": [output state vector]
			},
		...
		],
		"components": [Component Gates]
	},
	"Equiv Phase" : [Global Phase Ignorance Option]
}
```


Users need to insert content on part indicated with square brackets `[ ]`.
(for now,  field value of `"Spec Type"` and `"max_inst_allowed"` can be ignored). Specifically, required item are as follows:
- `"ID"` : (Any) name for the synthesis problem.	

- `"qreg_size"` : Specify number of qubits of circuit to synthesize in integer type value. This corresponds to $N$ in our definition of synthesis problem (see Section~3 in our paper for detail)


- `"spec_object"` :  Put list of input/output statevectors that circuit should satisfy. This corresponds to $E$ in our definition. Insert inputstate vector at `[input state vector]` and output statevector at `[output state vector]`. Valid form of statevector value should be written by sequence of (float) number in string value, with delimiter by comma `,`. For example, followings are valid form of in/output statevector (note that $0.7071..  = \frac{1}{\sqrt{2}}$).

	```
	"0,  0,  0,  0,  1, 0, 0, 0 " 

	" 0.70710678118,  0,  0,  0,  0, 0, 0, 0.70710678118 "
	```

	Each string value denotes $\ket{100},  \frac{1}{\sqrt{2}}(\ket{000} + \ket{111})$.
	In case of complex numbers, following Python's representation it should be written as `a+b j` where `a` and `b` is some float number (e.g, `0.25-0.25j`). 

- `"components"` : Specify component gates used for circuit presentation. This corresponds to $\mathcal{G}$ in our definition. At `[Component Gates]`, list name of component gates with delimiter `,` in string value. For example, to specify Hadamard gate and CNOT gate as component gate, insert following :
    ```
    "H, CNOT"
    ```
    Conceptually, our synthesis algorithm generally adopts any of quantum gate. However, efficient interface for user specifying component gate is  yet in development.
    
    Instead in an ad-hoc way, users may specify their component gate by adding items on dictionary variable of name `TO_CIRQ_GATE` at `/qsyn/set_synthesis.py` in format of
    \[
        \texttt{[name of gate]: [cirq.Gate instance] }.
    \]



- `"Equiv Phase"` : Whether to ignore global phase when validating circuit. Give in json boolean type value `true` to ignore global phase, otherwise `false`.

For example of this `.json` file template, refer any of file under `/benchmark` directory.

After generating such `.json` file, put it under `/benchmark` folder.
Then run following command to start synthesis:

```
python3 qsyn/run_single.py --benchmark [your file name]  --mode [Mode]
```
where `[your file name]` is name of `.json` file newly generated.



### Reproducing Qiskit's Transpile

Following command will show those transpiled circuits in printed text:
```
python3 transpile_qc.py
```

## Contact

Kang, Chan Gu (e-mail : changukang@korea.ac.kr)
