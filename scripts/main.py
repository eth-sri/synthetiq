import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from analyzer import main_analysis
import os
import json
from qiskit import QuantumCircuit, QuantumRegister
import subprocess
import warnings
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42; matplotlib.rcParams['ps.fonttype'] = 42

source_file = "data/times.csv"
# source_file = "data/times_paper.csv"
source_file_optimization = "data/optimization.csv"
source_output = 'data/output'
output_folder = "data/paper"
output_folder_qsyn = "qsyn/output"
mosca_data_path = 'data/baselines/mosca'
gheorghiu_data_path = 'data/baselines/gheorgiu'

kang_counts = {
    'bit_measure': 0,
    'cluster': 13,
    'draper': 5,
    'flip': 6,
    'GHZ_by_iSWAP': 8,
    'GHZ_by_QFT': 7,
    'GHZ_from_100': 6,
    'GHZ_Game': 8,
    'indexed_bell': 4,
    'M_valued': 7,
    'QFT': 7,
    'teleportation': 4,
    'three_superpose': 3,
    'toffoli_by_sqrt_X': 5,
    'W_four': 7,
    'W_orthog': 8,
    'W_phased': 7,
    '1': 0,
    '2': 0,
    '3': 0,
    '4': 0,
    '5': 0,
}

df = pd.read_csv(source_file, sep='\t', header=None)
df.columns = ['filename', 'time', 'std', 'n_runs', 'successful_runs', 'optimal_runs']

def add_rccx(circuit, c1, c2, t1):
    circuit.h(t1)
    circuit.cx(t1, c1)
    circuit.t(c1)
    circuit.cx(c2, t1)
    circuit.cx(c2, c1)
    circuit.t(c1)
    circuit.tdg(t1)
    circuit.cx(c2, t1)
    circuit.tdg(t1)
    circuit.cx(c2, c1)
    circuit.cx(c2, t1)
    circuit.cx(c1, t1)
    circuit.cx(c2, c1)
    circuit.cx(c2, t1)
    circuit.h(c1)
    circuit.cx(c1, t1)
    circuit.cx(t1, c1)
    return circuit

def add_ciswap(circuit, c1, c2, t1):
    circuit.cx(t1,c2)
    circuit.h(t1)
    circuit.cx(t1,c2)
    circuit.cx(c2,c1)
    circuit.t(c2)
    circuit.tdg(c1)
    circuit.cx(t1,c2)
    circuit.cx(c2,c1)
    circuit.tdg(t1)
    circuit.cx(t1,c2)
    circuit.t(c1)
    circuit.cx(c2,c1)
    circuit.cx(t1,c2)
    circuit.h(t1)
    circuit.cx(c2,c1)
    circuit.cx(t1,c2)
    circuit.cx(t1,c2)

def controlled_circuit(circuit, implementations):
    num_qubits = circuit.num_qubits

    new_qr = QuantumRegister(num_qubits + 1 + 1, name="q")

    controlled_circ = QuantumCircuit(new_qr, name="controlled_circuit")

    for inst, qargs, _ in circuit.data:
        gate = inst
        target_qubits = [new_qr[q.index + 1] for q in qargs]

        if gate.name in implementations:
            if gate.name == 't' or gate.name == "tdg":  # Replace 't' with the key you're using in the implementations dictionary for the T gate
                implementation_qubits = [new_qr[0]] + target_qubits + [new_qr[-1]]
            else:
                implementation_qubits = [new_qr[0]] + target_qubits

            implementation = implementations[gate.name]
            controlled_circ.compose(implementation, qubits=implementation_qubits, inplace=True)

        else:
            controlled_gate = gate.control(1)
            controlled_circ.append(controlled_gate, [new_qr[0]] + target_qubits)

    return controlled_circ

def save_circuit_to_qasm(circuit, filename):
    qasm_str = circuit.qasm()
    with open(filename, 'w') as file:
        file.write(qasm_str)

def table4():
    outputs = []
    os.makedirs(os.path.join(output_folder, '62'), exist_ok=True)
    for folder in os.listdir(os.path.join(source_output, '62')):
        output = dict()
        output['operator'] = folder
        _, _, _, _, best_t_depth_circ, _ = main_analysis(os.path.join(source_output, '62', folder))
        output['t_depth'] = best_t_depth_circ.t_depth
        output['qubits'] = best_t_depth_circ.circuit.num_qubits

        controlled_implementations = {
            "cx": QuantumCircuit.from_qasm_file("data/baselines/controlled/ccx.qasm"),
            "t": QuantumCircuit.from_qasm_file("data/baselines/controlled/ct.qasm"),
            "h": QuantumCircuit.from_qasm_file("data/baselines/controlled/ch.qasm"),
            "s": QuantumCircuit.from_qasm_file("data/baselines/controlled/cs.qasm"),
            "sdg": QuantumCircuit.from_qasm_file("data/baselines/controlled/cs_inverse.qasm"),
            "tdg": QuantumCircuit.from_qasm_file("data/baselines/controlled/ct_inverse.qasm"),
        }

        if folder == 'rcccx':
            output['baseline_t_depth'] = 8 # from paper https://arxiv.org/abs/1508.03273
        elif folder == 'cct':
            circuit = QuantumCircuit(5)
            add_rccx(circuit, 0, 1, 3)
            # add controlled t gate
            add_rccx(circuit, 2,3,4)
            circuit.t(4)
            add_rccx(circuit, 2,3,4)
            add_rccx(circuit, 0, 1, 3)
            output['baseline_t_depth'] = circuit.depth(lambda gate: gate[0].name in ['t', 'tdg'])
        elif folder == 'cciswap':
            circuit = QuantumCircuit(5) # using https://threeplusone.com/pubs/on_gates.pdf
            add_rccx(circuit, 0, 1, 4)
            add_ciswap(circuit, 2,3,4)
            add_rccx(circuit, 0, 1, 4)
            output['baseline_t_depth'] = circuit.depth(lambda gate: gate[0].name in ['t', 'tdg'])
        elif folder == 'csqrtswap':
            circuit = QuantumCircuit(2)
            circuit.cx(1,0)
            circuit.h(1)
            circuit.t(1)
            circuit.cx(0,1)
            circuit.tdg(1)
            circuit.tdg(0)
            circuit.h(1)
            circuit.cx(1,0)
            circuit.sdg(1)
            circuit.s(0)
            controlled_circ = controlled_circuit(circuit, controlled_implementations)
            save_circuit_to_qasm(controlled_circ, f"data/input/baseline.qasm")
            # running simplification for honest comparison
            output_subprocess = subprocess.check_output(f"bin/main_resynth baseline.qasm", shell=True)
            output_ = output_subprocess.decode("utf-8").split("T-depth: ")[1].split(" -> ")[1]
            # remove the baseline file
            os.remove(f"data/input/baseline.qasm")
            int_tdepth = int(output_)
            output['baseline_t_depth'] = int_tdepth
        else:
            circuit = QuantumCircuit(2)
            circuit.cx(1,0)
            circuit.h(1)
            circuit.tdg(1)
            circuit.cx(0,1)
            circuit.tdg(1)
            circuit.h(1)
            circuit.sdg(1)
            circuit.s(0)
            circuit.h(1)
            circuit.cx(1,0)
            circuit.sdg(0)
            controlled_circ = controlled_circuit(circuit, controlled_implementations)
            save_circuit_to_qasm(controlled_circ, f"data/input/baseline.qasm")
            # running simplification for honest comparison
            output_subprocess = subprocess.check_output(f"bin/main_resynth baseline.qasm", shell=True)
            output_ = output_subprocess.decode("utf-8").split("T-depth: ")[1].split(" -> ")[1]
            os.remove(f"data/input/baseline.qasm")
            int_tdepth = int(output_)
            output['baseline_t_depth'] = int_tdepth

        best_t_depth_circ.circuit.qasm(formatted=True, filename=os.path.join(output_folder, '62', f'{folder}.qasm'))
        outputs.append(output)
    
    outputs = pd.DataFrame(outputs)
    outputs.to_csv(os.path.join(output_folder, '62', 'table4.csv'), index=False)
        
def optimization():
    csv = pd.read_csv(source_file_optimization, sep="\t", header=None)
    csv.columns = ["name", "time", "std", "runs", "success", "real_success"]
    csv["qbs"] = csv["name"].apply(lambda x: int(x.split("/")[3]))
    csv["parameter"] = csv["name"].apply(lambda x: x.split("/")[4])
    csv["value"] = csv["name"].apply(lambda x: x.split("/")[5] if not x.split("/")[5].endswith("txt") else 0)

    csv["value"] = csv["value"].astype(float)
    csv["matrix"] = csv["name"].apply(lambda x: x.split("/")[6] if len(x.split("/")) == 8 else x.split("/")[5])
    csv["type"] = csv["matrix"].apply(lambda x: x.split("_")[1])

    should_add = csv["real_success"] < 100
    numbers = np.maximum(np.array(csv["real_success"]), 1)
    csv["time"] = csv["time"] + should_add * ((600 - csv["time"] * csv["real_success"]) / numbers)

    best_values = csv.groupby(["parameter", "value"]).mean(numeric_only=True)["time"].reset_index()
    best_times = best_values.groupby("parameter").min()["time"].reset_index()
    best_values = best_values.merge(best_times, on=["parameter", "time"], how="inner")

    all_csv_parameters = None
    for parameter in best_values["parameter"].unique():
        csv_parameter = csv[csv["parameter"] == parameter].copy()
        csv_for_merging = csv_parameter[['matrix', 'qbs', 'value', "time"]].copy()
        # merge the time of the best value on each matrix
        parameter_value = best_values[best_values["parameter"] == parameter]["value"].values[0]
        csv_parameter = csv_parameter.merge(csv_for_merging[csv_for_merging["value"] == parameter_value], on=["matrix", "qbs"], how="inner")
        if all_csv_parameters is None:
            all_csv_parameters = csv_parameter
        else:
            all_csv_parameters = pd.concat([all_csv_parameters, csv_parameter])

    all_csv_parameters["relative"] = all_csv_parameters["time_x"] / all_csv_parameters["time_y"]

    grouped = all_csv_parameters[["relative", "parameter", "value_x"]].groupby(["parameter", "value_x"]).mean()
    output = dict()
    # get value optimal for each parameter
    for parameter in grouped.index.get_level_values(0).unique():
        # set parameter equal to the value that has the lowest relative time
        output[parameter] = grouped.loc[parameter].idxmin().values[0]

    naming_conversion = {
        'n-norm': r'n_{\text{norm}}',
        'pid': r'P_{\text{Id}}',
        'start-temp': 'T(0)',
        'n-start-gates': r'\ell',
        'iterations-factor': r'n_{\text{step}} / (n\ell)',
    }

    output = {naming_conversion[k]: v for k, v in output.items() if k in naming_conversion.keys()}
    
    os.makedirs(os.path.join(output_folder, '61'), exist_ok=True)
    json.dump(output, open(os.path.join(output_folder, "61", "parameters.json"), "w"), indent=4)

    csv_only_table2 = csv[csv["parameter"].isin(["", "no-resynth", "simple", "no-perms"])] 
    csv_normal = csv_only_table2[csv_only_table2["parameter"] == ""]
    csv_normal = csv_normal[["matrix", "qbs", "time"]]
    csv_only_merged = csv_only_table2.merge(csv_normal, on=["matrix", "qbs"], how="inner")
    csv_only_merged["relative"] = csv_only_merged["time_x"] / csv_only_merged["time_y"] 

    csv_only_merged = csv_only_merged[csv_only_merged["parameter"] != ""]
    grouped = csv_only_merged[["relative", "parameter", "type"]].groupby(["parameter", "type"]).mean()
    grouped.to_csv(os.path.join(output_folder, '61', "table2.csv"))

def fig3():
    synthetiq = df.copy()
    synthetiq = synthetiq[synthetiq["filename"].str.startswith("data/output/63/theirs/")]
    # synthetiq = synthetiq[synthetiq["filename"].str.startswith("output/62/theirs/")]
    synthetiq = synthetiq[np.logical_not(synthetiq["filename"].str.contains("count"))]
    synthetiq['benchmark'] = synthetiq["filename"].apply(lambda x: x.split("/")[4])
    # synthetiq['benchmark'] = synthetiq["filename"].apply(lambda x: x.split("/")[3])
    should_add = synthetiq["optimal_runs"] < 100
    synthetiq["time"] = synthetiq["time"] + should_add * (10800 / synthetiq["optimal_runs"] - synthetiq["time"])

    kang_data = []
    for benchmark in synthetiq['benchmark']:
        with open(os.path.join(output_folder_qsyn, benchmark + '.txt')) as f:
            raw_output = f.read()
            failure = ''
            if 'Elapsed Time :' in raw_output:
                time = float(raw_output.split('Elapsed Time :')[1].split('\n')[0])
                num_gates = kang_counts[benchmark]
            elif 'TIMEOUT' in raw_output:
                failure = 'T/O'
                time = 24 * 3600
                num_gates = 0
            else:
                failure = 'N/A'
                time = 24 * 3600
                num_gates = 0
            kang_data.append({
                'benchmark': benchmark,
                'time': time,
                'num_gates': num_gates,
                'failure': failure
            })
    kang = pd.DataFrame(kang_data)
    
    for benchmark in synthetiq['filename']:
        benchmark_name = benchmark.replace("data/output/", "")
        # benchmark_name = benchmark.replace("output/", "").replace('62/', '63/')
        min_num_gates = 100000000
        for file in os.listdir(os.path.join(source_output, benchmark_name)):
            # get the amount of non-empty lines in the file
            with open(os.path.join(source_output, benchmark_name, file)) as f:
                num_gates = sum(1 for line in f if line.strip()) - 3 # -3 because of the header
            min_num_gates = min(min_num_gates, num_gates)
        synthetiq.loc[synthetiq['filename'] == benchmark, 'num_gates'] = min_num_gates
    benchmarks = synthetiq['benchmark'].unique()
    # change the order of the benchmarks list to follow the order of (kang['num_gates'], synthetiq['num_gates'])
    zero_to_high = lambda x: 1000 if x == 0 else x
    benchmarks = sorted(benchmarks, key=lambda x: (zero_to_high(kang[kang['benchmark'] == x]['num_gates'].values[0]), 
                                                   synthetiq[synthetiq['benchmark'] == x]['num_gates'].values[0]))
    # sort the dataframes by the benchmarks list
    synthetiq = synthetiq.set_index('benchmark').loc[benchmarks].reset_index()
    kang = kang.set_index('benchmark').loc[benchmarks].reset_index()

    renames = {
            '1': 'CX',
            '2': 'sqrt_SWAP',
            '3': 'Permutation',
            '4': 'sqrt_iSWAP',
            '5': 'Oracle',
    }
    synthetiq['benchmark'] = synthetiq['benchmark'].apply(lambda x: renames.get(x, x))
    kang['benchmark'] = kang['benchmark'].apply(lambda x: renames.get(x, x))
    speedups = np.array(kang["time"]) / np.array(synthetiq["time"])
    speedups = [speedup for i, speedup in enumerate(speedups) if kang["failure"][i] == '']
    fig, ax = plt.subplots(figsize=(12, 5), ncols=2, width_ratios=[5, 1])
    current_palette = sns.color_palette("Dark2", 2)
    sns.boxplot(speedups, width=0.8, color=current_palette[0], ax=ax[1])
    ax[1].set_yscale('log')
    ax[1].set_xticklabels([])
    ax[1].set_xticks([])
    ax[1].set_title("Speedup")
    # make the lines surrounding the plot disappear
    sns.despine(left=True, bottom=True)
    #make the background grey
    ax[1].set_facecolor((0.95, 0.95, 0.95))
    # make the title to the left
    ax[1].title.set(x=0.05)
    # background color very light grey
    ax[0].set_facecolor((0.95, 0.95, 0.95))
    ax[0].title.set(x=0)
    ax[0].set_title("Gate Count")
    sns.despine(left=True, bottom=True)

    data = pd.DataFrame({'Source': ['Synthetiq' for i in range(len(synthetiq))] + ['Kang et al.' for i in range(len(kang))],
                        'count': np.concatenate([synthetiq['num_gates'], kang['num_gates']]), 
                        'label': np.concatenate([benchmarks, benchmarks])}
    )

    sns.set_palette(current_palette)
    sns.barplot(data=data, y="count", x="label", hue="Source", edgecolor=(0.95, 0.95, 0.95), linewidth=0, ax=ax[0])

    # rotate the xticks 45 degrees, making sure that they still appear under the barplots at the bottom
    ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation=45, ha='right', fontsize=12)

    ax[0].set_ylabel("")
    ax[0].set_xlabel("")
    ax[0].set_xlim(-0.5, 21.5)
    ax[0].xaxis.labelpad = 20
    ax[0].text(16.05, 2.7, "T/O", rotation=90, fontsize=11, color=current_palette[1], weight="bold")
    ax[0].text(17.05, 4.2, "N/A", rotation=90, fontsize=11, color=current_palette[1], weight="bold")
    ax[0].text(18.05, 4.2, "N/A", rotation=90, fontsize=11, color=current_palette[1], weight="bold")
    ax[0].text(19.05, 6.2, "T/O", rotation=90, fontsize=11, color=current_palette[1], weight="bold")
    ax[0].text(20.05, 6.7, "N/A", rotation=90, fontsize=11, color=current_palette[1], weight="bold")
    ax[0].text(21.05, 9.2, "T/O", rotation=90, fontsize=11, color=current_palette[1], weight="bold")
    # set title font size
    ax[0].title.set_fontsize(16)
    ax[0].tick_params(axis='y', labelsize=16)
    # set the legend font size
    ax[0].legend(fontsize=16)
    ax[1].title.set_fontsize(16)

    # set xticks font size
    ax[1].tick_params(axis='y', labelsize=16)

    ax[1].set_ylim(10 ** -1.2, 1.5 * 10 ** 5)

    plt.tight_layout() 

    fig.savefig('data/paper/fig_3.pdf', bbox_inches='tight')


def table5():
    results = []
    os.makedirs(os.path.join(output_folder, '63', 'stackexchange'), exist_ok=True)
    for folder in os.listdir(os.path.join(source_output, '63', 'stackexchange')):
        output = dict()
        full_folder = os.path.join(source_output, '63', 'stackexchange', folder)
        # search in df for the time
        element = df[df['filename'] == full_folder + '/'].to_dict('records')[0]
        output['name'] = {
            '1': 'CX',
            '2': '\sqrt{SWAP}',
            '3': 'Permutation',
            '4': '\sqrt{iSWAP}',
            '5': 'Oracle',
        }[folder]
        output['time'] = element['time']
        _, _, _, circ, _, _ = main_analysis(full_folder)
        output['num_qubits'] = circ.circuit.num_qubits
        output['expensive_count'] = circ.count
        circ.circuit.qasm(formatted=True, filename=os.path.join(output_folder, '63', 'stackexchange', f'{folder}.qasm'))
        results.append(output)

    output = pd.DataFrame(results)
    output.to_csv(os.path.join(output_folder, '63', 'table5.csv'), index=False)

def parse_mosca(filename):
    with open(filename) as f:
        raw_output = f.read()
        if 'Path = ' not in raw_output:
            return None, None
        time = raw_output.split('Execution time = ')[1].split('\n')[0]
        time = float(time)
        num_gates = len(raw_output.split('Path = ')[1].split('\n')[0].split(':'))

    return time, num_gates

def parse_gheorgiu(filename):
    with open(filename) as f:
        raw_output = f.read()
        if 'OUT:' not in raw_output:
            return None, None
        time = raw_output.split('\n')[-1].split('Took: ')[1].split(' ms')[0]
        time = float(time) / 1000
        num_gates = int(raw_output.split('T-depth: ')[1].split('\n')[0])

    return time, num_gates


def table6():
    results = []
    os.makedirs(os.path.join(output_folder, '64', 'comparison'), exist_ok=True)
    for folder in os.listdir(os.path.join(source_output, '64', 'comparison', 'tcount')):
        output = dict()
        full_folder = os.path.join(source_output, '64', 'comparison', 'tcount', folder)
        # search in df for the time
        try:
            element = df[df['filename'] == full_folder + '/'].to_dict('records')[0]
        except IndexError:
            warnings.warn(f"Could not find an optimal implementation for {folder} in the dataframe. This is caused due to variation between the runs. Running the command for this again should fix this.")
            continue
        if element['optimal_runs'] < 100: # need to add to time until we are at 12 hours
            element['time'] = 12 * 3600 / element['optimal_runs']
        element['time'] = round(element['time'], 2)
        
        output['name'] = folder
        output['tcount_time'] = element['time']
        _, _, _, best_t_count_circ, _, _ = main_analysis(full_folder)
        output['tcount'] = best_t_count_circ.t_count
        best_t_count_circ.circuit.qasm(formatted=True, filename=os.path.join(output_folder, '64', 'comparison', f'{folder}_tcount.qasm'))
        time_baseline, gates_baseline = parse_mosca(os.path.join(mosca_data_path, 'comparison', f'{folder}.txt'))
        output['baseline_time_count'] = time_baseline
        output['baseline_gates_count'] = gates_baseline
        if time_baseline is not None:
            output['tcount_speedup'] = time_baseline / element['time']
        
        full_folder = os.path.join(source_output, '64', 'comparison', 'tdepth', folder)
        # search in df for the time
        try:
            element = df[df['filename'] == full_folder + '/'].to_dict('records')[0]
        except IndexError:
            warnings.warn(f"Could not find an optimal implementation for {folder} in the dataframe. This is caused due to variation between the runs. Running the command for this again should fix this.")
            continue
        if element['optimal_runs'] < 100: # need to add to time until we are at 12 hours
            element['time'] = 12 * 3600 / element['optimal_runs']
        element['time'] = round(element['time'], 2)
        
        output['name'] = folder
        output['tdepth_time'] = element['time']
        _, _, _, _, best_t_depth_circ, _ = main_analysis(full_folder)
        output['tdepth'] = best_t_depth_circ.t_depth

        time_baseline, gates_baseline = parse_gheorgiu(os.path.join(gheorghiu_data_path, 'comparison', f'{folder}.txt'))
        output['baseline_time_depth'] = time_baseline
        output['baseline_gates_depth'] = gates_baseline
        if time_baseline is not None:
            output['tdepth_speedup'] = time_baseline / element['time']
        
        best_t_depth_circ.circuit.qasm(formatted=True, filename=os.path.join(output_folder, '64', 'comparison', f'{folder}_tdepth.qasm'))

        results.append(output)
    
    results = pd.DataFrame(results)
    results.to_csv(os.path.join(output_folder, '64', 'comparison', 'table6.csv'), index=False)

def load_data_fig4():
    csv = df.copy()
    csv.columns = ["folder", "time", "std", "runs", "succ_runs", "real_success"]
    # csv = csv[csv["folder"].str.startswith("output/63/permutations/")]
    csv = csv[csv["folder"].str.startswith("data/output/64/permutations/")]
    csv.drop(["std", "runs", "succ_runs"], axis=1, inplace=True)
    csv["kind_theirs"] = csv["folder"].apply(lambda x: x.split("/")[-4])
    csv["kind"] = csv["folder"].apply(lambda x: x.split("/")[-3])
    csv["theirs"] = csv["kind_theirs"].isin(["tcount", "tdepth"])
    csv["kind"] = np.where(csv["theirs"], csv["kind_theirs"], csv["kind"])
    csv["number"] = csv["folder"].apply(lambda x: int(x.split("/")[-2]))
    csv.drop(["folder"], inplace=True, axis=1)
    # adjust the time where number equals to set to the one where csv["theirs"] is true
    csv = csv.merge(csv[csv["theirs"]][["time", "number", "kind"]], how="left", on=["number", "kind"])
    csv["time_ours"] = csv["time_x"]
    csv["time"] = np.where(csv["time_y"].isna(), csv["time_x"], csv["time_y"])
    csv.drop(["time_x", "time_y"], inplace=True, axis=1)
    csv = csv[np.logical_not(csv["theirs"])]
    
    groups_same = {
        0: [0],
        1: [1,2,5,6,13,22,29],
        2: [3,4,7,8,11,12,15,16,19,20,23,24,27,28],
        3: [9,10,14,17,18,21,25,26]
    }

    # invert the groups, make each elemtn in the list a single key
    groups = {v: k for k, vs in groups_same.items() for v in vs}

    normal_data = {
        0: (0, 0), #T-count, T-depth
        1: (7, 3),
        2: (8, 3),
        3: (15, 7)
    }
    csv["tcount"] = csv["number"].apply(lambda x: normal_data[groups[x]][0])
    csv["tdepth"] = csv["number"].apply(lambda x: normal_data[groups[x]][1])
    csv.drop(["real_success"], axis=1, inplace=True)
    
    mosca_data = []
    for i in range(0, 30):
        with open(os.path.join(mosca_data_path, 'permutations', str(i) + ".txt")) as f:
            tcount = None
            time = 0
            for line in f.readlines():
                if line.startswith("Path :  "):
                    tcount = len(line.split(","))
                if line.startswith("Execution time :  "):
                    time = float(line.strip()[18:])
            mosca_data.append({
                "number": i,
                "time": time,
                "tcount": tcount,
                "kind": "tcount"
            })
            
    mosca = pd.DataFrame(mosca_data)
    
    gheorgiou_data = [{
                "number": 0,
                "time": None,
                "tdepth": None,
                "kind": "tdepth"
            }]
    file = os.path.join(gheorghiu_data_path, 'permutations.txt')
    with open(file) as f:
        current_number = None
        for i, line in enumerate(f.readlines()):
            if i % 4 == 0:
                current_number = int(line.strip())
            elif i % 4 == 2:
                if "T-depth" in line:
                    tdepth = int(line.split(" ")[1])
                else:
                    tdepth = None
            elif i % 4 == 3:
                time = float(line.split(" ")[1]) / 1000
                gheorgiou_data.append({
                    "number": current_number,
                    "time": time,
                    "tdepth": tdepth,
                    "kind": "tdepth"
                })

    gheorgiou = pd.DataFrame(gheorgiou_data)
    
    csv["origin"] = "ours"
    gheorgiou["origin"] = "gheorgiou"
    mosca["origin"] = "mosca"
    # merge all three
    all_data = pd.concat([csv, gheorgiou, mosca])
    
    csv_tcount = csv[csv["kind"] == "tcount"]
    csv_tcount_copy = csv_tcount.copy()
    csv_tcount_copy.sort_values(["tcount", "time"], inplace=True)
    csv_tcount_copy.reset_index(inplace=True)
    other_numbering = {
        el["number"]: i for i, el in enumerate(csv_tcount_copy.to_dict(orient="records"))
    }
    
    csv_tdepth = csv[csv["kind"] == "tdepth"]
    csv_tdepth_copy = csv_tdepth.copy()
    csv_tdepth_copy.sort_values(["tdepth", "time"], inplace=True)
    csv_tdepth_copy.reset_index(inplace=True)
    other_numbering_tdepth = {
        el["number"]: i for i, el in enumerate(csv_tdepth_copy.to_dict(orient="records"))
    }
    
    tcount_data = all_data[all_data["kind"] == "tcount"].copy()
    tdepth_data = all_data[all_data["kind"] == "tdepth"].copy()
    tcount_data["number"] = tcount_data["number"].apply(lambda x: other_numbering[x])
    tdepth_data["number"] = tdepth_data["number"].apply(lambda x: other_numbering_tdepth[x])
    
    return tcount_data, tdepth_data

def fig4():
    tcount, tdepth = load_data_fig4()

    mosca = tcount[tcount["origin"] == "mosca"].copy()
    synthetiq = tcount[tcount["origin"] == "ours"].copy()

    # order them by number
    mosca.sort_values(["number"], inplace=True)
    synthetiq.sort_values(["number"], inplace=True)

    # reorder them further such that mosca tcount goes from 0 to max to nan
    mosca.sort_values(["tcount"], inplace=True)
    # now make synthetiq get the same order
    synthetiq["number"] = synthetiq["number"].apply(lambda x: mosca[mosca["number"] == x]["number"].values[0])
    synthetiq.sort_values(["number"], inplace=True)

    fig, ax = plt.subplots(figsize=(14, 3), width_ratios=[12, 2], ncols=2)
    if any(synthetiq['time'] == 0):
        warnings.warn(f"There are {np.count_nonzero(synthetiq['time'] == 0)} permutations for which no optimal implementation was found within the time limit. This is caused due to the variation between runs. Running the command again should fix this.")
    speedups = np.array(mosca["time"]) / np.array(synthetiq["time"])
    # plt.scatter([1 for i in range(len(speedups))], speedups, label='Synthetiq', color='blue')
    current_palette = sns.color_palette("Dark2", 2)
    sns.boxplot([speedup for speedup in speedups if speedup > 0], width=0.8, color=current_palette[0], ax=ax[1])
    ax[1].set_yscale('log')
    ax[1].set_xticklabels([])
    ax[1].set_xticks([])
    ax[1].set_title("Speedup")
    # make the lines surrounding the plot disappear
    sns.despine(left=True, bottom=True)
    #make the background grey
    ax[1].set_facecolor((0.95, 0.95, 0.95))
    # make the title to the left
    ax[1].title.set(x=0.05)
    ax[1].set_title("Speedup")

    mosca_counts = [count if count is not None else 0 for count in mosca["tcount"]]
    synthetiq_counts = list(synthetiq["tcount"])

    data = pd.DataFrame({'Source': ['Synthetiq' for i in range(len(synthetiq_counts))] + ['Mosca et al.' for i in range(len(mosca_counts))],
                        'count': synthetiq_counts + mosca_counts, 
                        'label': [i for i in range(len(synthetiq_counts))] + [i for i in range(len(mosca_counts))]}
    )

    ax[0].set_title("T-Count")
    ax[0].title.set(x=0)

    ax[0].set_facecolor((0.95, 0.95, 0.95))

    sns.set_palette(current_palette)
    sns.barplot(data=data, y="count", x="label", hue="Source", edgecolor=(0.95, 0.95, 0.95), linewidth=0, ax=ax[0])

    ax[0].set_xticks([])
    ax[0].set_ylabel("")
    ax[0].set_xlabel("")
    ax[0].set_xlim(-0.5, 29.5)
    ax[0].xaxis.labelpad = 20
    plt.tight_layout()
    for i in range(22, 30):
        ax[0].text(i, 7, "T/O", rotation=90, fontsize=10, color=current_palette[1], weight="bold")

    ax[0].text(-0.3, 0.2, "0", rotation=0, fontsize=10, color=current_palette[0], weight="bold")  

    # set title font size
    ax[0].title.set_fontsize(16)
    ax[0].tick_params(axis='y', labelsize=16)
    ax[0].tick_params(axis='x', labelsize=16)
    # set the legend font size
    ax[0].legend(fontsize=16)
    ax[1].title.set_fontsize(16)

    # set xticks font size
    ax[1].tick_params(axis='y', labelsize=16)

    fig.tight_layout()

    ax[1].set_ylim(10 ** -1.2, 1.5 * 10 ** 5)

    fig.savefig('data/paper/fig_4a.pdf', bbox_inches='tight')

    synthetiq = tdepth[tdepth["origin"] == "ours"].copy()
    gheorgiou = tdepth[tdepth["origin"] == "gheorgiou"].copy()

    # order them by number
    gheorgiou.sort_values(["number"], inplace=True)
    synthetiq.sort_values(["number"], inplace=True)

    if any(synthetiq['time'] == 0):
        warnings.warn(f"There are {np.count_nonzero(synthetiq['time'] == 0)} permutations for which no optimal implementation was found within the time limit. This is caused due to the variation between runs. Running the command again should fix this.")

    speedups = np.array(gheorgiou["time"]) / np.array(synthetiq["time"])

    gheorgiou["ours_tdepth"] = synthetiq["tdepth"].tolist()
    gheorgiou.sort_values(["ours_tdepth", "tdepth"], inplace=True)

    fig, ax = plt.subplots(figsize=(14, 3), ncols=2, width_ratios=[12, 2])
    # plt.scatter([1 for i in range(len(speedups))], speedups, label='Synthetiq', color='blue')
    current_palette = sns.color_palette("Dark2", 2)
    sns.boxplot([speedup for speedup in speedups if speedup > 0], width=0.8, color=current_palette[0], ax=ax[1])
    ax[1].set_yscale('log')
    ax[1].set_xticklabels([])
    ax[1].set_xticks([])
    ax[1].set_title("Speedup")
    # make the lines surrounding the plot disappear
    sns.despine(left=True, bottom=True)
    #make the background grey
    ax[1].set_facecolor((0.95, 0.95, 0.95))
    # make the title to the left
    ax[1].title.set(x=0.05)

    ax[0].set_facecolor((0.95, 0.95, 0.95))
    ax[0].title.set(x=0)
    ax[0].set_title("T-Depth")
    sns.despine(left=True, bottom=True)

    gheorgiou_counts = [count if count is not None else 0 for count in gheorgiou["tdepth"]]
    synthetiq_counts = list(gheorgiou["ours_tdepth"])

    data = pd.DataFrame({'Source': ['Synthetiq' for i in range(len(synthetiq_counts))] + ['Gheorgiou et al.' for i in range(len(gheorgiou_counts))],
                        'count': synthetiq_counts + gheorgiou_counts, 
                        'label': [i for i in range(len(synthetiq_counts))] + [i for i in range(len(gheorgiou_counts))]}
    )

    sns.set_palette(current_palette)
    sns.barplot(data=data, y="count", x="label", hue="Source", edgecolor=(0.95, 0.95, 0.95), linewidth=0, ax=ax[0])

    ax[0].set_xticks([])
    ax[0].set_ylabel("")
    ax[0].set_xlabel("")
    ax[0].set_xlim(-0.5, 29.5)
    ax[0].xaxis.labelpad = 20
    plt.tight_layout()
    for i in range(23, 30):
        ax[0].text(i, 3.3, "T/O", rotation=90, fontsize=10, color=current_palette[1], weight="bold")
        
    for i in range(17, 22):
        ax[0].text(i, 1.3, "T/O", rotation=90, fontsize=10, color=current_palette[1], weight="bold")
        
    # ax[0].text(7, 1.3, "T/O", rotation=90, fontsize=10, color=current_palette[1], weight="bold")

    ax[0].text(0.1, 0.2, "T/O", rotation=90, fontsize=10, color=current_palette[1], weight="bold")
    ax[0].text(-0.2, 0.3, "0", rotation=0, fontsize=10, color=current_palette[0], weight="bold")


    # set title font size
    ax[0].title.set_fontsize(16)
    ax[0].tick_params(axis='y', labelsize=16)
    ax[0].tick_params(axis='x', labelsize=16)
    # set the legend font size
    ax[0].legend(fontsize=14)
    ax[1].title.set_fontsize(16)

    # set xticks font size
    ax[1].tick_params(axis='y', labelsize=16)
    ax[1].set_ylim(10 ** -1.2, 1.5 * 10 ** 5)

    fig.tight_layout()

    fig.savefig('data/paper/fig_4b.pdf', bbox_inches='tight')

def extract_results_baseline_epsilon(file):
    new_lines = []
    with open(file, "r") as f:
        current_time = 0
        start_time = False
        for line in f.readlines():
            if ".txt" in line or line.startswith("#### T-count:"):
                if start_time:
                    new_lines.append(f"Took: {current_time}\n")
                new_lines.append(line)
                start_time = False
                current_time = 0
            elif line.startswith("Took: "):
                current_time += int(line.split(" ")[1])
                start_time = True

    results = []
    for i in range(0, len(new_lines), 3):
        name = new_lines[i].split('_')[0]
        if name == 'rzk':
            name = 'rz'
        if name == 'crzk':
            name = 'crz'
        k_value = new_lines[i].split('_')[1]
        epsilon = float(new_lines[i].split('_')[2][:-5])
        ms = int(new_lines[i+1][5:-1])
        try:
            t_count = int(new_lines[i+2][13:-6])
        except ValueError:
            t_count = None

        results.append((name + '_' + k_value, epsilon, ms / 1000, t_count))

    dataframe = pd.DataFrame(results, columns=["name", "epsilon", "time", "t_count"])
    return dataframe

def table7():
    results = []
    baseline_results = extract_results_baseline_epsilon(os.path.join(gheorghiu_data_path, 'epsilon.txt'))
    os.makedirs(os.path.join(output_folder, '65'), exist_ok=True)
    for folder in os.listdir(os.path.join(source_output, '65')):
        for subfolder in os.listdir(os.path.join(source_output, '65', folder)):
            output = dict()
            full_folder = os.path.join(source_output, '65', folder, subfolder)
            # search in df for the time
            element = df[df['filename'] == full_folder + '/'].to_dict('records')[0]
            _, t_count, _, _, best_t_depth_circ, _ = main_analysis(full_folder)
            output['name'] = subfolder
            output['epsilon'] = {'005': 0.05, '10_2': 0.01, '10_3': 0.001}[folder]
            baseline_element = baseline_results[(baseline_results['name'] == subfolder) & (baseline_results['epsilon'] == output['epsilon'])]
            if len(baseline_element) == 0:
                output['baseline_time'] = None
                output['baseline_tcount'] = None
            else:
                output['baseline_time'] = baseline_element['time'].values[0]
                output['baseline_tcount'] = baseline_element['t_count'].values[0]
            if len(t_count) == 0:
                output['time'] = None
                output['tcount'] = None
            else:
                min_t_count = min(t_count)
                number_min_t_count = len(t_count[t_count == min_t_count])
                if element['successful_runs'] >= 100:
                    output['time'] = float(element['time'] * element['successful_runs'] / number_min_t_count)
                else:
                    output['time'] = float(3600 / number_min_t_count)
                output['tcount'] = min_t_count
                if output['baseline_time'] is not None and output['time'] is not None:
                    output['speedup'] = output['baseline_time'] / output['time']

                best_t_depth_circ.circuit.qasm(formatted=True, filename=os.path.join(output_folder, '65', f'{subfolder}_{str(output["epsilon"])}.qasm'))

            results.append(output)
    results = pd.DataFrame(results)
    results.to_csv(os.path.join(output_folder, '65', 'table7.csv'), index=False)

def ccx_tdepth():
    circ = QuantumCircuit(3)
    circ.h(2)
    circ.t(0)
    circ.t(1)
    circ.t(2)
    circ.cx(1, 0)
    circ.cx(2, 1)
    circ.cx(0, 2)
    circ.tdg(1)
    circ.cx(0, 1)
    circ.tdg(0)
    circ.tdg(1)
    circ.t(2)
    circ.cx(2, 1)
    circ.cx(0, 2)
    circ.h(2)
    circ.cx(1, 0)
    return circ.to_gate()

def ccx_tcount():
    circ = QuantumCircuit(3)
    circ.h(2)
    circ.cx(1, 2)
    circ.tdg(2)
    circ.cx(0, 2)
    circ.t(2)
    circ.cx(1, 2)
    circ.t(1)
    circ.tdg(2)
    circ.cx(0, 2)
    circ.cx(0, 1)
    circ.t(2)
    circ.h(2)
    circ.t(0)
    circ.tdg(1)
    circ.cx(0, 1)
    return circ.to_gate()

def cirq_carry_gate(ccxg):
    # https://github.com/quantumlib/Cirq/blob/main/examples/basic_arithmetic.py
    carry_circ = QuantumCircuit(4)
    carry_circ.append(ccxg, [1, 2, 3])
    carry_circ.cx(1, 2)
    carry_circ.append(ccxg, [0, 2, 3])
    return carry_circ.to_gate()

def cirq_uncarry_gate(ccxg):
    # https://github.com/quantumlib/Cirq/blob/main/examples/basic_arithmetic.py
    carry_circ = QuantumCircuit(4)
    carry_circ.append(ccxg, [0, 2, 3])
    carry_circ.cx(1, 2)
    carry_circ.append(ccxg, [1, 2, 3])
    return carry_circ.to_gate()

def build_adder(n, carry_gate, uncarry_gate):
    # https://github.com/quantumlib/Cirq/blob/main/examples/basic_arithmetic.py
    def sum_gate():
        sum_circ = QuantumCircuit(3)
        (c0, a, b) = (sum_circ.qubits[0], sum_circ.qubits[1], sum_circ.qubits[2])
        sum_circ.cx(a, b)
        sum_circ.cx(c0, b)
        return sum_circ.to_gate()

    a = QuantumRegister(n, name = "a")
    b = QuantumRegister(n, name = "b")
    c = QuantumRegister(n, name = "c")

    circuit = QuantumCircuit(a, b, c)

    for i in range(n-1):
        circuit.append(carry_gate, [c[i], a[i], b[i], c[i+1]])
    circuit.append(sum_gate(), [c[n-1], a[n-1], b[n-1]])
    for i in range(n-2, -1, -1):
        circuit.append(uncarry_gate, [c[i], a[i], b[i], c[i+1]])
        circuit.append(sum_gate(), [c[i], a[i], b[i]])

    return circuit

def plot_fig6(size, nb_gates_cirq, nb_gates_ours, label):
    fig, ax = plt.subplots(figsize=(3,3))
    current_palette = sns.color_palette("Dark2", 2)
    sns.despine(fig, left=True, bottom=True)
    ax.set_facecolor((0.95, 0.95, 0.95))
    ax.plot(size, nb_gates_cirq, color=current_palette[1])
    ax.plot(size, nb_gates_ours, color=current_palette[0])
    ax.title.set(x=0.05)
    ax.set_title(label)
    ax.set_xlabel('# Qubits', fontsize=12)
    fig.tight_layout()
    os.makedirs(os.path.join(output_folder, 'fig6'), exist_ok=True)
    fig.savefig(os.path.join(output_folder, 'fig6', f'{label}.pdf'))

def fig6():
    _, _, _, _, optimal_carry_t, _ = main_analysis(os.path.join(source_output, '66', 'carry_relative'))
    _, _, _, _, _, optimal_carry_cx = main_analysis(os.path.join(source_output, '66', 'carry_relative_cnot'))
    optimal_uncarry_t = optimal_carry_t.circuit.copy().inverse().to_gate()
    optimal_carry_t = optimal_carry_t.circuit.to_gate()
    optimal_uncarry_cx = optimal_carry_cx.circuit.copy().inverse().to_gate()
    optimal_carry_cx = optimal_carry_cx.circuit.to_gate()

    x = np.arange(2, 21, 1)
    y_tdepth = []
    for nq in x:
        adder_circ = build_adder(nq, cirq_carry_gate(ccx_tdepth()), cirq_uncarry_gate(ccx_tdepth()))
        decomp_circ = adder_circ.decompose().decompose()
        #print(decomp_circ)
        depth = decomp_circ.depth(lambda g: g[0].name in ['t', 'tdg'])
        y_tdepth.append(depth)
    
    y_tcount = []
    for nq in x:
        adder_circ = build_adder(nq, cirq_carry_gate(ccx_tcount()), cirq_uncarry_gate(ccx_tcount()))
        decomp_circ = adder_circ.decompose().decompose()
        #print(decomp_circ)
        count_ops = decomp_circ.count_ops()
        count = count_ops['t'] + count_ops['tdg']
        y_tcount.append(count)

    y_cx_depth = []
    for nq in x:
        adder_circ = build_adder(nq, cirq_carry_gate(ccx_tcount()), cirq_uncarry_gate(ccx_tcount()))
        decomp_circ = adder_circ.decompose().decompose()
        #print(decomp_circ)
        depth = decomp_circ.depth(lambda g: g[0].name in ['cx'])
        y_cx_depth.append(depth)
    
    y_cx_count = []
    for nq in x:
        adder_circ = build_adder(nq, cirq_carry_gate(ccx_tcount()), cirq_uncarry_gate(ccx_tcount()))
        decomp_circ = adder_circ.decompose().decompose()
        #print(decomp_circ)
        count_ops = decomp_circ.count_ops()
        count = count_ops['cx']
        y_cx_count.append(count)

    y_tdepth_optimal = []
    y_tcount_optimal = []
    for nq in x:
        adder_circ = build_adder(nq, optimal_carry_t, optimal_uncarry_t)
        decomp_circ = adder_circ.decompose()
        #print(decomp_circ)
        depth = decomp_circ.depth(lambda g: g[0].name in ['t', 'tdg'])
        count_ops = decomp_circ.count_ops()
        count = count_ops['t'] + count_ops['tdg']
        y_tdepth_optimal.append(depth)
        y_tcount_optimal.append(count)

    y_cx_depth_optimal = []
    y_cx_count_optimal = []

    for nq in x:
        adder_circ = build_adder(nq, optimal_carry_cx, optimal_uncarry_cx)
        decomp_circ = adder_circ.decompose()
        #print(decomp_circ)
        depth = decomp_circ.depth(lambda g: g[0].name in ['cx'])
        count_ops = decomp_circ.count_ops()
        count = count_ops['cx']
        y_cx_depth_optimal.append(depth)
        y_cx_count_optimal.append(count)

    plot_fig6(x, y_tdepth, y_tdepth_optimal, 'T-Depth')
    plot_fig6(x, y_tcount, y_tcount_optimal, 'T-Count')
    plot_fig6(x, y_cx_depth, y_cx_depth_optimal, 'CX-Depth')
    plot_fig6(x, y_cx_count, y_cx_count_optimal, 'CX-Count')

    
if __name__ == "__main__":
    optimization()
    fig3()
    table4()
    table5()
    table6()
    fig4()
    table7()
    fig6()