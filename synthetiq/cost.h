#ifndef DEF_COST
#define DEF_COST

#include "circuit.h"
#include "randomhelper.h"
#include "partialMatrix.h"
#include <Eigen/Dense>

typedef Eigen::Array<bool, Eigen::Dynamic, 1> BoolArray;

class EqualityComputer {
   public:
        double normalizedEqualityCost(GateCircuit& circ, QubitIndependentPartialMatrix& matrix_obj, CircuitHelper& ch);
        virtual double normalizedEqualityCost(GateCircuit& circ, PartialMatrix& matrix_obj, CircuitHelper& ch); 
        virtual std::shared_ptr<EqualityComputer> clone() = 0;
};

class ExactEqualityComputer : public EqualityComputer {
    public:
        ExactEqualityComputer() {};
        ExactEqualityComputer(ExactEqualityComputer const& other) : tolerance(other.tolerance) {};
        std::shared_ptr<EqualityComputer> clone();
        ExactEqualityComputer(double tolerance): tolerance(tolerance) {};
        double normalizedEqualityCost(GateCircuit& circ, QubitIndependentPartialMatrix& matrix_obj, CircuitHelper& ch);
        double normalizedEqualityCost(GateCircuit& circ, PartialMatrix& matrix_obj, CircuitHelper& ch); 
    private:
        double tolerance = 1e-6;
};
 
class FroebeniusCostComputer: public EqualityComputer {
   public:
        FroebeniusCostComputer();
        FroebeniusCostComputer(FroebeniusCostComputer const& other) {};
        std::shared_ptr<EqualityComputer> clone();
        double normalizedEqualityCost(GateCircuit& circ, PartialMatrix& matrix_obj, CircuitHelper& ch); 
};

class SimpleFroebeniusCostComputer: public EqualityComputer {
   public:
        SimpleFroebeniusCostComputer();
        SimpleFroebeniusCostComputer(SimpleFroebeniusCostComputer const& other) {};
        std::shared_ptr<EqualityComputer> clone();
        double normalizedEqualityCost(GateCircuit& circ, PartialMatrix& matrix_obj, CircuitHelper& ch); 
};

class PerformanceComputer{
    public:
        virtual double normalizedPerformanceCost(GateCircuit& c, CircuitHelper& ch) = 0;
        virtual std::shared_ptr<PerformanceComputer> clone() = 0;
};

class GatesSumComputer: public PerformanceComputer{
    public:
        GatesSumComputer(GatesSumComputer const& other);
        GatesSumComputer();
        std::shared_ptr<PerformanceComputer> clone(); 
        virtual double normalizedPerformanceCost(GateCircuit& c, CircuitHelper& ch);
};

#endif