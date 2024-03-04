#ifndef DEF_CIRCUIT_HELPER
#define DEF_CIRCUIT_HELPER
#include <vector>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <iostream>
#include <string>
#include <fstream>
#include <filesystem>
#include <regex>
#include <limits>
#include "randomhelper.h"
#include "gate.h"

class CircuitHelper {
    public:
        CircuitHelper(int nb_qbs = 1, std::string basic_gate_folder="data/gates/CliffordT", std::string composite_gate_folder="data/gates/composite_gates", std::string read_gate_folder="data/gates/read_gates");
        std::vector<std::shared_ptr<Gate>> all_gates;
        std::vector<std::shared_ptr<Gate>> readable_gates;
        std::vector<std::shared_ptr<Gate>> basic_gates;
        std::vector<std::shared_ptr<Gate>> composite_gates;
        std::vector<std::vector<std::shared_ptr<Gate>>> basic_gates_by_name;
        std::vector<std::vector<std::shared_ptr<Gate>>> composite_gates_by_name;
        std::shared_ptr<Gate> id_gate;
        std::shared_ptr<Gate> invertGate(std::shared_ptr<Gate> gate);
        std::string basic_gate_folder;
        std::string composite_gate_folder;
        std::string read_gate_folder;
        int nb_qbs;
        double max_cost_basic;
        double max_cost_all;
    private:
        bool isAlreadyPresent(std::string name, std::vector<int> acting_qubits);
        void readBasicGateFolder(std::string folder);
        void readCompositeGateFolder(std::string folder, bool read_folder=false);
};

#endif
