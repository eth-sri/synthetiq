import cirq
if __name__ == "qsyn.module_restrain":
    from qsyn.util.utils import *
else:
    from util.utils import *

def restrain_to_sp_modules(ss ,  m : cirq.Circuit, targ_diff : tuple)-> bool:
    for g_op in m.all_operations():
        if g_op.gate in ss.component_prior[NC] or g_op.gate in ss.component_prior[ENTANGLE]: return True
    return False

def restrain_to_phase_modules(ss, m:cirq.Circuit, targ_diff :tuple )-> bool:
    if includes_phase_others_bool(ss, m) :
        return True
    if superpose_open_and_closed(ss, m):
        return True
    return False

def restrain_to_bool_modules(ss , curr_state, m : cirq.Circuit, targ_diff : tuple)-> bool:
    bool_li = list()
    if targ_diff[1]:
        for q in [ss.working_qubit[i] for i in targ_diff[1]]:
            for i in range(len(m)):
                temp = m.operation_at(q,i)
                if temp: bool_li.append(True)
                else: bool_li.append(False)
        if not any(bool_li): 
            return False
    if curr_state.is_curr_out_arr_all_classic and is_oddly_sp(ss,m):
        return False
    if is_not_all_phase(ss, m):
        if is_all_bool(ss, m) or superpose_open_and_closed(ss, m):
            return True
        return False
    else :
        return False

def includes_phase_others_bool(ss, m):
    at_least_one_phase = False
    for g in m.all_operations():
        if g.gate in ss.component_prior[PHASING]:
            at_least_one_phase = True
        elif g.gate not in ss.component_prior[BOOL]:
            return False
    return at_least_one_phase 

def is_oddly_sp(ss, m : cirq.Circuit):
    cnt = 0 
    for g in m.all_operations():
        if g.gate in ss.component_prior[NC] +  ss.component_prior[ENTANGLE] :
            cnt+=1
    return cnt == 1 or cnt%2==1  ## HARD CODED

def is_all_bool(ss, m : cirq.Circuit):
    
    if len(ss.component_prior[BOOL]) == 0:
        return False

    for x in m.all_operations():
        if x.gate not in ss.component_prior[BOOL]:return False
    return True

def is_not_all_phase(ss, m: cirq.Circuit):
    if len(ss.component_prior[PHASING]) == 0:
        return True

    for x in m.all_operations():
        if x.gate not in ss.component_prior[PHASING]:return True
    return False

def superpose_open_and_closed(ss, m:cirq.Circuit):
    if len(m) == 1:
        return False
    first_moment = m[0]
    last_moment  = m[-1]
    mapping_pair = dict()

    for first_m_gate_op in first_moment:
        if  first_m_gate_op.gate not in (ss.component_prior[NC] + ss.component_prior[ENTANGLE]):
            return False
    for last_m_gate_op in last_moment:
        if  last_m_gate_op.gate not in (ss.component_prior[NC] + ss.component_prior[ENTANGLE]):
            return False
            
    for first_m_gate_op in first_moment :

        for qubit_of_first_m in first_m_gate_op.qubits:
            last_m_gate_op = last_moment.operation_at(qubit_of_first_m)
            if last_m_gate_op != None:
                break

        if first_m_gate_op.gate in (ss.component_prior[NC] + ss.component_prior[ENTANGLE]):
            if last_m_gate_op == None :
                return False
            elif isinstance(first_m_gate_op, cirq.ControlledOperation) and isinstance(last_m_gate_op, cirq.ControlledOperation): # such as CH,...
                if set(first_m_gate_op.sub_operation.qubits) == set(last_m_gate_op.sub_operation.qubits) :
                    continue
                else : return False
            elif (last_m_gate_op.qubits) ==  (first_m_gate_op.qubits) and last_m_gate_op.gate in (ss.component_prior[NC] + ss.component_prior[ENTANGLE]):
                continue
            else : 
                return False
        elif last_m_gate_op :
            # if last_m_gate_op.gate not in ss.component_prior[BOOL]:
            return False
    return True

def constrained_module_evolve(ss, m: cirq.Circuit):
    # because generating module of level 4 is so expensive,
    # we generated the only that will be only needed, from level 3
    # m : module of level 3 = num of gate is 3
    # print(m)
    assert superpose_open_and_closed(ss, m) or is_all_bool(ss,m)
    classic_gate_op = list()
    phasing_gate_op = list()
    for gate_op in ss.gate_operations:
        if gate_op.gate in ss.component_prior[BOOL]:
            classic_gate_op.append(gate_op)
        if gate_op.gate in ss.component_prior[PHASING]:
            phasing_gate_op.append(gate_op)
    to_return = list()
    if is_all_bool(ss,m):
        for c_gate_op in classic_gate_op :
            copied_m = m.copy()
            copied_m.append(c_gate_op)
            to_return.append(copied_m)
    elif superpose_open_and_closed(ss, m):
        for phase_gate_op in phasing_gate_op :
            copied_m = m.copy()
            copied_m.insert(-1, phase_gate_op)
            if superpose_open_and_closed(ss, copied_m):
                # print(copied_m)
                # print("=======")
                # input()
                to_return.append(copied_m)
    return to_return
