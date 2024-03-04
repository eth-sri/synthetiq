import numpy as np

from state_search import *
from set_synthesis import *
from guide import *
from util.utils import *
from state import *
from bfs_baseline import exectue_synthesis
# from qsyn.state_search import *
# from qsyn.main import *
# from qsyn.guide import *
# from qsyn.util.utils import *
# from qsyn.state import *
from timeit import default_timer as timer




if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', type = str, default="Ours", help="Put name of synthesis algorithm." 
                                                                    + "\nIt can be one of \'Ours\', \'Ours_no_prune\', \'Base\', \'Base_no_prune\'."
                                                                    + "\nFor detaield specification of each algorithm mode, see our paper.")
    args = parser.parse_args()
    benchmarks = [
        'three_superpose'
        ,'M_valued'
        ,'GHZ_from_100'
        ,'GHZ_by_iSWAP'
        ,'GHZ_by_QFT'
        ,'GHZ_Game'
        ,'W_orthog'
        ,'W_phased'
        ,'W_four'
        ,'cluster'
        ,'bit_measure'
        ,'flip'
        ,'teleportation'
        ,'indexed_bell'
        ,'toffoli_by_sqrt_X'
        ,'QFT'
        ,'draper'
    ]

    res_time_collector = dict()
    if args.mode == "Ours":
        for idx, benchmark in enumerate(benchmarks):
            spec = load_benchmark(benchmark)
            print(spec)
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
                                    timeout = 3600)
            current_time = timer()
            print("Elapsed Time", current_time - start_time)
            print("========================")            
            if res_qc :
                print(f"Number of Gate In synthesized Circuit : {count_gate(res_qc)}")   
                if current_time - start_time >= 3600 :
                    res_time_collector[idx] = [benchmark, "TIME OUT", "NA", "NA"]
                else:
                    res_time_collector[idx] = [benchmark, current_time - start_time, count_gate(res_qc), len(res_modules)]
            else :
                if current_time - start_time >= 3600 :
                    res_time_collector[idx] = [benchmark, "TIME OUT", "NA", "NA"]
                else :
                    res_time_collector[idx] = [benchmark, "NA", "NA", "NA"]
            print("Synthesized")
            print(res_qc)       
            
    elif args.mode == "Ours_no_prune":
        for idx, benchmark in enumerate(benchmarks):
            spec = load_benchmark(benchmark)
            print(spec)
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
                                    timeout = 3600)
            current_time = timer()
            print("Elapsed Time", current_time - start_time)
            print("========================")            
            if res_qc :
                print(f"Number of Gate In synthesized Circuit : {count_gate(res_qc)}")   
                if current_time - start_time >= 3600 :
                    res_time_collector[idx] = [benchmark, "TIME OUT"]
                else:
                    res_time_collector[idx] = [benchmark, current_time - start_time]
            else :
                if current_time - start_time >= 3600 :
                    res_time_collector[idx] = [benchmark, "TIME OUT"]
                else :
                    res_time_collector[idx] = [benchmark, "NA"]
            print("Synthesized")
            print(res_qc)
    elif args.mode == "Base":
        for idx, benchmark in enumerate(benchmarks):
            elapsed_time, synth_circuit , _ = exectue_synthesis(benchmark, 3600, True, False)
            res_time_collector[idx] = [benchmark, elapsed_time]
    elif args.mode == "Base_no_prune":
        for idx, benchmark in enumerate(benchmarks):
            elapsed_time, synth_circuit , _ = exectue_synthesis(benchmark, 3600, False, False)
            res_time_collector[idx] = [benchmark, elapsed_time]
    else : 
        print("invalid paramter for --mode. Only put one of Ours / Ours_no_prune / Base / Base_no_prune")
        exit()
    
    if args.mode == "Ours":
        print("=================================")       
        print("Experiment Result : Ours")
        print("=================================")     
        print("{:<8} {:<20} {:<15} {:<5} {:<5}".format('Number','Benchmark','Time(s)', "#O", "#M"))
        for k, v in res_time_collector.items():
            label, syn_time, num_gate, num_module = v
            print("{:<8} {:<20} {:<15} {:<5} {:<5}".format(k+1, label , syn_time if isinstance(syn_time, str) else round(syn_time,2) , num_gate, num_module   ))
    else :
        print("=================================")        
        print(f"Experiment Result : {args.mode}")
        print("================================")
        print("{:<8} {:<20} {:<15} ".format('Number','Benchmark','Time(s)'))
        for k, v in res_time_collector.items():
            label, syn_time = v
            print("{:<8} {:<20} {:<15}".format(k+1, label , syn_time if isinstance(syn_time, str) else round(syn_time,2) ))
    