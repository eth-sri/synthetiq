import numpy as np

from state_search import *
from set_synthesis import *
from guide import *
from util.utils import *
from state import *
from bfs_baseline import exectue_synthesis
from timeit import default_timer as timer
import numpy as np


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--benchmark', type=str,required=True, help="Put ID of the benchmark." 
                                                                    + "\nThe synthesis problem should be presented in \'.json\' file in proper format, under \'/benchmark\' folder."
                                                                    + "\nFor the format of the file, see instruction in artifact manual or github readme.md.")
    parser.add_argument('--mode', type = str, default="Ours", help="Put name of synthesis algorithm." 
                                                                    + "\nIt can be one of \'Ours\', \'Ours_no_prune\', \'Base\', \'Base_no_prune\'."
                                                                    + "\nFor detaield specification of each algorithm mode, see our paper.")
    parser.add_argument("--timeout", type=int, default=3600, help="Timeout for synthesis in seconds. Default is 3600s (1 hour).")
    

    args = parser.parse_args()
    spec = load_benchmark(args.benchmark)
    if args.mode == "Ours" :
        start_time = timer() 
        module_gate_num = spec.qreg_size
        ss = StateSearchSynthesis(spec=spec,
                                prepare_naive_module_gen = True,
                                entangle_template = True,
                                module_gate_num = module_gate_num)
        res_qc, res_modules = ss.modular_search(init_qc=None, 
                                        rule_depth=2, 
                                        concrete_criterion=False,
                                        naive_module_gen=True,
                                        entangle_template = True,
                                        no_pruning = False,
                                        no_mutate =  True,
                                        start_time=start_time,
                                        timeout = args.timeout)
        current_time = timer()
        elapsed_time = current_time - start_time

    elif args.mode == "Ours_no_prune":  
        start_time = timer() 
        module_gate_num = spec.qreg_size
        ss = StateSearchSynthesis(spec=spec,
                                prepare_naive_module_gen = True,
                                entangle_template = True,
                                module_gate_num = module_gate_num)
        res_qc, res_modules = ss.modular_search(init_qc=None, 
                                        rule_depth=2, 
                                        concrete_criterion=False,
                                        naive_module_gen=True,
                                        entangle_template = True,
                                        no_pruning = True,
                                        no_mutate =  True,
                                        start_time=start_time,
                                        timeout = args.timeout)
        current_time = timer()
        elapsed_time = current_time - start_time
    elif args.mode == "Base":
        elapsed_time, res_qc , _ = exectue_synthesis(args.benchmark, args.timeout, True, False)
    elif args.mode == "Base_no_prune":
        elapsed_time, res_qc , _ = exectue_synthesis(args.benchmark, args.timeout, False, False)
    else : 
        print("invalid paramter for --mode. Only put one of Ours / Ours_no_prune / Base / Base_no_prune")
        exit()

    print("================================")
    print("Synthesis Result")
    print(f"Benchmark : {args.benchmark}")
    print(f"Mode : {args.mode}")
    print("================================")

    if res_qc == None:
        print("None Found")
    else :
        print("Synthesized QC")
        print(res_qc)
        if args.mode == "Ours":
            print("Stacked Moduels")
            print(res_modules)
        print(f"Elapsed Time :  {round(elapsed_time,2)}")
