#ifndef DEF_CIRCUIT
#define DEF_CIRCUIT
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
#include "circuithelper.h"
#include "matrix_computer.h"

class GateCircuit {
    public:
        GateCircuit(CircuitHelper& ch, MatrixComputerType matrix_computer_type=Binary);
        GateCircuit(int nb_gates, int nb_q, CircuitHelper& ch, MatrixComputerType matrix_computer_type=Binary);
        GateCircuit(GateCircuit const& other);
        GateCircuit(std::vector<std::shared_ptr<Gate>> gates, int nb_q, CircuitHelper& ch, MatrixComputerType matrix_computer_type=Binary);
        ~GateCircuit() {}
        int nbElements(); //must return total possible number of gates (including Id)(if representation changes and we don't have a fixed number of gates anymore, beware of PerfCost computation, may not be normalized anymore)
        int performanceNormalizationCst();
        void readGateFromQasmInputLine(std::string line);
        void readFromInput(std::string filename); //from cin
        void readFromInputNonQuasm(std::string filename, int nb_l_to_ignore, int nb_qbs_in_cirs); //from cin
        void readFromIfStream(std::ifstream& file);
        std::string print_qasm();
        void calculateCost();
        void updateCost(std::shared_ptr<Gate> new_gate, std::shared_ptr<Gate> old_gate);
        double getCost();
        int getNbNonIdGates();
        Eigen::MatrixXcd toMatrix();
        const std::vector<std::shared_ptr<Gate>> getGates();
        bool mutate(RandomHelper& rh, double proportional_prob=1.0, double proba_id=0.2, double proba_name=1.0);
        bool mutate(RandomHelper& rh, int position, double proportional_prob=1.0, double proba_id=0.2, double proba_name=1.0);

        void undoMutation();
        void applyMutationAgain(GateCircuit& other);
        void placeGateAt(int position, std::shared_ptr<Gate> gate); //replaces gate at pos with new gate, assumes qbs of size 2
        void changeQubits(std::vector<int> qbs_order);
        void stopMatrixComputer();
        void startMatrixComputer();
        void expandCompositeGates();
        int getCount(std::vector<std::string> gate_names);
        int getDepth(std::vector<std::string> gate_names);
        void invert();
        void rotate();
        void initializeMatrixComputer();
        std::shared_ptr<CircuitHelper> getCircuitHelper();
        std::shared_ptr<GateCircuit> clone();

        int nb_qbs;
        std::shared_ptr<Gate> old_gate;
        std::shared_ptr<Gate> new_gate;
        int pos_mutation;

    private:
        std::vector<std::shared_ptr<Gate>> list_gates;
        std::shared_ptr<MatrixComputer> matrixComputer;
        std::shared_ptr<CircuitHelper> ch;
        int nb_non_id_gates = 0;
        double cost = 0;
        MatrixComputerType matrix_computer_type;
        bool calculating_matrix_computer = true;
};

#endif
