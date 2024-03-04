#ifndef DEF_RESYNTH
#define DEF_RESYNTH
#include "circuit.h"


class Resynthesize {
    public:
        Resynthesize(int max_gate_mult=12, bool watch_depth=false, std::vector<std::string> depth_gates = {"t", "tdg"}) : max_gate_mult(max_gate_mult), watch_depth(watch_depth), depth_gates(depth_gates) {};
        void run(std::shared_ptr<GateCircuit> circuit, CircuitHelper& ch, bool second_time=false);
    private:
        int change(std::shared_ptr<GateCircuit> circuit, int start_index, CircuitHelper& ch, bool second_time=false);
        bool commutes(std::shared_ptr<Gate> gate1, std::shared_ptr<Gate> gate2);
        int NbBeforeCommutes(int index, std::shared_ptr<GateCircuit> circuit);
        bool priorityChange(int index, std::shared_ptr<GateCircuit> circuit, bool second_time=false);
        int max_gate_mult;
        bool watch_depth;
        std::vector<std::string> depth_gates;
};

#endif