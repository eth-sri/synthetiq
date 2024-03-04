import cirq
import logging
from typing import List, Tuple, Union
import argparse
if __name__ == "bfs_baseline":
    from synthesis_spec.specification import Spec
    from synthesis_spec.utils import IOpairs
    from synthesis import Synthesis, BFSSynthesis
    from set_synthesis import load_benchmark
else : 
    from synthesis_spec.specification import Spec
    from synthesis_spec.utils import IOpairs
    from synthesis import Synthesis, BFSSynthesis
    from set_synthesis import load_benchmark
    from util.utils import count_gate
    
from timeit import default_timer as timer

def exectue_synthesis(benchmark_id: str, timeout: int,  involution_prune: bool, identity_prune: bool) -> Tuple[Union[str, int], cirq.Circuit]:
    loaded_spec = load_benchmark(benchmark_id)
    # logging.info("Printed Spec :" + "\n" + str(loaded_spec))
    syn_app = BFSSynthesis(spec=loaded_spec)
    # logging.info("Component_Prior: " + "\n" +
                #  pprint.pformat(syn_app.component_prior, indent=4))
    # logging.info("(Spanned) Gate Operatios Set" +
                #  "\n" + str(syn_app.spanned_set))
    # logging.info(f"# Of Spanned Gate Operation : {len(syn_app.spanned_set)}")
    # logging.info(f"The size of Search Space is {syn_app.get_search_space_size()}")
    # <=> number of nodes in the search tree except for root.

    print(loaded_spec)
    elapsed_time, synth_circuit = syn_app.prune_bfs_enumeration_synthesis(timeout=timeout, involution_prune=involution_prune, identity_prune=identity_prune)
    return elapsed_time, synth_circuit, loaded_spec


if __name__ == '__main__':
    # prepare arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--benchmark', type=str,required=True, help="ID of the benchmark")
    parser.add_argument('--timeout', type=int, required=True)
    parser.add_argument('--involution-prune', default=False, action='store_true')
    parser.add_argument('--identity-prune', default=False, action='store_true')
    parser.add_argument('--save_in_json', default=False, action='store_true')
    args = parser.parse_args()
    
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(filename="[LOG]"+str(args.benchmark),
                        filemode='a',
                        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.INFO)
    logging.info("Start Main Procedure")
    elapsed_time, synth_circuit, loaded_spec = exectue_synthesis(args.benchmark, args.timeout, args.involution_prune, args.identity_prune)
    print(f"Elapsed Time :  {elapsed_time}")
    print("========================")
    if synth_circuit and args.save_in_json:
        json_string = cirq.to_json(synth_circuit)
        print('JSON string:')
        print(json_string)
    if synth_circuit:
        print(f"Number of Gate In synthesized Circuit : {count_gate(synth_circuit)}")
    print("Synthesized")
    print(synth_circuit)

    logging.info(f"\n{str(synth_circuit)}")
