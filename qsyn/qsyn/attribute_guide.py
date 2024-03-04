import numpy as np

from state_search import *
from set_synthesis import *
from guide import *
from util.utils import *
from state import *

from timeit import default_timer as timer



if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--benchmark', type=str, required=True, help="ID of the benchmark")
    parser.add_argument('--concrete_criterion', default=False, action='store_true')
    parser.add_argument('--naive_module_gen', default=False, action='store_true')
    parser.add_argument('--entangle_template', default=False, action='store_true')
    parser.add_argument('--main_module_gen', default=False, action='store_true')
    parser.add_argument('--no_pruning', default=False, action='store_true')
    parser.add_argument('--no_mutate', default=False, action='store_true')
    parser.add_argument('--save_in_json', default=False, action='store_true')
    parser.add_argument('--module_gate_num', type = int, default=None)
    parser.add_argument('--timeout', type = int, default=300)

    args = parser.parse_args()
    assert not (args.entangle_template == True and args.naive_module_gen == False)




    spec = load_benchmark(args.benchmark)
    print(spec)
    start_time = timer() 

    if args.module_gate_num == None :
        module_gate_num = spec.qreg_size
    else :
        module_gate_num  = args.module_gate_num
    if args.main_module_gen == True:
        ss = StateSearchSynthesis(spec=spec,
                                prepare_naive_module_gen = True,
                                entangle_template = True,
                                module_gate_num = module_gate_num)
    else:
        ss = StateSearchSynthesis(spec=spec,
                                prepare_naive_module_gen = args.naive_module_gen,
                                entangle_template = args.entangle_template,
                                module_gate_num = module_gate_num)
    if args.main_module_gen == True:
        res_qc, res_modules = ss.modular_search(init_qc=None, 
                                        rule_depth=2, 
                                        concrete_criterion=args.concrete_criterion,
                                        naive_module_gen=True,
                                        entangle_template = True,
                                        no_pruning = args.no_pruning,
                                        no_mutate =  args.no_mutate,
                                        start_time=start_time,
                                        timeout = args.timeout)
    else:
        res_qc, res_modules = ss.modular_search(init_qc=None, 
                                        rule_depth=2, 
                                        concrete_criterion=args.concrete_criterion,
                                        naive_module_gen=args.naive_module_gen,
                                        entangle_template = args.entangle_template,
                                        no_pruning = args.no_pruning,
                                        no_mutate =  args.no_mutate,
                                        start_time=start_time,
                                        timeout = args.timeout)

    current_time = timer()
    print("Elapsed Time", current_time - start_time)
    print("========================")

    if res_qc and args.save_in_json:
        json_string = cirq.to_json(res_qc)
        print('JSON string:')
        print(json_string)
    if res_qc :
        print(f"Number of Gate In synthesized Circuit : {count_gate(res_qc)}")   
    print("Synthesized")
    print(res_qc)


