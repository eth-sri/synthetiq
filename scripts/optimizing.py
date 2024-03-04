import os
import numpy as np
import subprocess
import time
import argparse

parser = argparse.ArgumentParser(description='Run optimization')
parser.add_argument('--cores', type=int, default=64, help='Number of cores to use')

args = parser.parse_args()

store_file = "optimization"
input_folder = lambda qbs: f"data/input/61/{qbs}qbs"
normal_input_folder = lambda qbs: f"61/{qbs}qbs"
get_gate_count = lambda file: int(file.split("_")[0])
get_tcount = lambda file: int(file.split("_")[-3])

def get_command(input_file, qbs, parameter, parameter_value, tcount=None, extra_args=""):
    file = f"{store_file}/{qbs}/{parameter}/{parameter_value}/{input_file.split('/')[-1]}"
    file = file.replace("//", "/")
    command = f"./bin/main {input_file} -o {file} --{parameter} {parameter_value} --threads 1 --circuits 100 --time 600 --save --times-file data/optimization.csv --no-save-circuits"
    if tcount is not None:
        command += f" -tc {tcount}"
    else:
        command += " -gc 1000"
    command += " " + extra_args
    if parameter != "beta":
        command += " --beta 0"
    return command

all_commands = []
for qbs in [2, 3, 4]:
    parameters = {
        "start-temp": np.array([0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.16, 0.2]),
        "iterations-factor": np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]),
        "n-norm": np.array([30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]),
        "n-start-gates-factor": np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]),
        "pid": np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5]),
        "simple": [""],
        "no-perms": [""],
        "no-resynth": [""],
        "": [""]
    }

    for file in os.listdir(input_folder(qbs)):
        for parameter in parameters:
            for parameter_value in parameters[parameter]:                    
                n_gates = get_gate_count(file)
                tcount = get_tcount(file)
                if parameter == "n-start-gates-factor":
                    n_start_gates = int(n_gates * parameter_value)
                    command = get_command(f"{normal_input_folder(qbs)}/{file}", qbs, "n-start-gates", n_start_gates, tcount)
                elif parameter == "pid":
                    n_start_gates = int(n_gates * 3)
                    command = get_command(f"{normal_input_folder(qbs)}/{file}", qbs, parameter, parameter_value, tcount, "--n-start-gates " + str(n_start_gates))
                elif parameter == "pcomp":
                    command = get_command(f"{normal_input_folder(qbs)}/{file}", qbs, parameter, parameter_value, extra_args="-cg composite_rccx")
                else:
                    command = get_command(f"{normal_input_folder(qbs)}/{file}", qbs, parameter, parameter_value, tcount=tcount)
                all_commands.append(command)

cores = np.arange(0,args.cores)

running_tasks = [None for _ in cores]
current_task = 0

while current_task < len(all_commands):
    for i, task in enumerate(running_tasks):
        if task is None:
            running_tasks[i] = subprocess.Popen(f"taskset -c {i} " + all_commands[current_task], shell=True)
            current_task += 1
        elif task.poll() is not None:
            running_tasks[i] = None

    time.sleep(1)