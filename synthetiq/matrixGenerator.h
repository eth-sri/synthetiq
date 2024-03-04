#ifndef DEF_GENERATOR
#define DEF_GENERATOR

#include "partialMatrix.h"
#include <Eigen/Dense>
#include <Eigen/Core>

class MatrixGenerator {
    public:
        MatrixGenerator() {};
        PartialMatrix addAncilla(PartialMatrix& matrix, int n_ancilla, std::string name="matrix");
        PartialMatrix addDirtyQubits(PartialMatrix& matrix, int n_dirty, std::string name="matrix");
};

#endif