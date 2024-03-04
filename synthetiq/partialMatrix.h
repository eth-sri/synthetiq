#ifndef DEF_PARTIAL
#define DEF_PARTIAL
#include <vector>
#include <Eigen/Dense>
#include "circuit.h"

typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> BoolMatrix;

class BaseConstraint {
    public:
        BaseConstraint();
        BaseConstraint(int n_qubits) : n_qubits(n_qubits) {};
        virtual int nConstraints() = 0;
        int getNQubits();
        int n_qubits;
};



class PartialMatrix : public BaseConstraint {
    public:
        PartialMatrix();
        PartialMatrix(std::string filename);
        PartialMatrix(PartialMatrix const& other);
        PartialMatrix(GateCircuit& circuit, CircuitHelper& ch, std::string name="matrix");
        PartialMatrix(Eigen::MatrixXcd matrix, BoolMatrix cover, std::string name="matrix");
        PartialMatrix clone();
        int nConstraints();
        void write(std::string filename);
        Eigen::MatrixXcd getMatrix();
        Eigen::ArrayXXcd getArray();
        BoolMatrix getCover();
        void calculateInitialization();

        std::string name;
        double squared_norm;
        Eigen::MatrixXcd matrix;
        Eigen::ArrayXXcd array;
        BoolMatrix cover;        
};

class QubitIndependentPartialMatrix {
    public:
        QubitIndependentPartialMatrix (PartialMatrix& matrix, bool use_independent_qbs=true, bool use_inverse=false);
        QubitIndependentPartialMatrix() {};
        QubitIndependentPartialMatrix(std::string filename, bool use_independent_qbs=true, bool use_inverse=false);
        QubitIndependentPartialMatrix(QubitIndependentPartialMatrix const& other);
        QubitIndependentPartialMatrix(GateCircuit& circuit, CircuitHelper& ch, std::string name="matrix", bool use_independent_qbs=true, bool use_inverse=false);
        QubitIndependentPartialMatrix(Eigen::MatrixXcd matrix, BoolMatrix cover, std::string name="matrix", bool use_independent_qbs=true, bool use_inverse=false);
        void calculateQubitIndependent();
        int getNQubits();
        bool addMatrix(PartialMatrix& matrix);
        std::shared_ptr<QubitIndependentPartialMatrix> clone();

        PartialMatrix original;
        std::vector<std::shared_ptr<PartialMatrix>> matrices;
        std::vector<std::vector<int>> qbs_info;
        std::vector<bool> inverse_info;
        bool use_inverse;
        bool use_independent_qbs;
};

#endif