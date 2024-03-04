#ifndef DEF_GATE
#define DEF_GATE
#include <vector>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <iostream>
#include <string>
#include <fstream>
#include <filesystem>
#include <regex>
#include <limits>

class Gate {
    public:
        Gate();
        Gate(std::string name, Eigen::MatrixXcd matrix, std::vector<int> acting_qubits, double cost);
        std::string name;
        Eigen::MatrixXcd matrix;
        int nb_qbs;
        double cost;
        virtual bool isBasicGate() = 0;
        std::vector<int> acting_qubits; // qubits on which the gate acts, e.g. H on 0th qb, or on 1st qb, ...
        std::string print();
        virtual std::vector<std::shared_ptr<Gate>> decomposeInBasicGates() = 0;
        bool equals(Gate& other_gate);
        static std::shared_ptr<Gate> findCorrectGate(std::vector<std::shared_ptr<Gate>> allowed_gates, std::string gate_name, std::vector<int> acting_qbs_gate);
};

class BasicGate : public Gate {
    public:
        BasicGate(std::string filename);
        BasicGate(std::string name, Eigen::MatrixXcd matrix, std::vector<int> acting_qubits, double cost);
        bool isBasicGate();
        std::vector<std::shared_ptr<Gate>> decomposeInBasicGates();
};

class CompositeGate : public Gate {
    public:
        CompositeGate(std::string filename, std::vector<std::shared_ptr<Gate>> allowed_gates);
        CompositeGate(std::string name, Eigen::MatrixXcd matrix, std::vector<int> acting_qubits, double cost, std::vector<std::shared_ptr<Gate>> decomposition);
        std::vector<std::shared_ptr<Gate>> decomposition;
        std::vector<std::shared_ptr<Gate>> decomposeInBasicGates();
        void readFromFileTxT(std::string filename, std::vector<std::shared_ptr<Gate>> allowed_gates);
        void readFromFileQasm(std::string filename, std::vector<std::shared_ptr<Gate>> allowed_gates);
        bool isBasicGate();
        void changeQubitsDecomposition(std::vector<int> qbs_change, std::vector<std::shared_ptr<Gate>> allowed_gates);
};

#endif