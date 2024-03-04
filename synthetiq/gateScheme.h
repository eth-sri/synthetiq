#ifndef DEF_GATE_SCHEME
#define DEF_GATE_SCHEME
#include <vector>
#include <Eigen/Dense>
#include "randomhelper.h"

class GateScheme {
    public:
        GateScheme(int min_start_gates=30, int max_start_gates=120, double beta=0.05, double min_factor=2.5, double max_factor=3.5);
        void reset();
        int getStartGates(RandomHelper& rh);
        void update(int gate_count_found);
        void setMinStartGates(int min_start_gates);
        void setMaxStartGates(int max_start_gates);
        void setStartBestGates(int start_best_gates);
        void setBeta(double beta);
        double getBeta();
        void setMinFactor(double min_factor);
        void setMaxFactor(double max_factor);
        int getMinStartGates();
        int getMaxStartGates();
        double getMinFactor();
        double getMaxFactor();
        int getMinGates();
        int getMaxGates();

    private:
        int min_start_gates;
        int max_start_gates;
        int min_gates;
        int max_gates;
        int current_best_gates;
        double beta;
        double min_factor;
        double max_factor;
        int start_best_gates;
};

#endif