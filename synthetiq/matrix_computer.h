#ifndef DEF_MATRIX_COMPUTER
#define DEF_MATRIX_COMPUTER
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

enum MatrixComputerType {Linear, Chunk, Binary};

class MatrixComputer {
    public:
        virtual void updateMatrix(int position, std::vector<std::shared_ptr<Gate>> list_gates) = 0;
        virtual void  calculateMatrix(std::vector<std::shared_ptr<Gate>> list_gates) = 0;
        virtual Eigen::MatrixXcd getMatrix() = 0;
};

class LinearMatrixComputer : public MatrixComputer {
    public:
        LinearMatrixComputer();
        void updateMatrix(int position, std::vector<std::shared_ptr<Gate>> list_gates);
        void  calculateMatrix(std::vector<std::shared_ptr<Gate>> list_gates);
        Eigen::MatrixXcd getMatrix();
        Eigen::MatrixXcd matrix;
};

class ChunkMatrixComputer : public MatrixComputer {
    public:
        ChunkMatrixComputer();
        ChunkMatrixComputer(int nb_gates, int nb_qbs);
        void updateMatrix(int position, std::vector<std::shared_ptr<Gate>> list_gates);
        void  calculateMatrix(std::vector<std::shared_ptr<Gate>> list_gates);
        void initializeChunks(int nb_gates, int nb_qbs);
        Eigen::MatrixXcd getMatrix();
        std::vector<Eigen::MatrixXcd> chunks;
        Eigen::MatrixXcd matrix;
        bool matrix_computed = false;
        int nb_gates;
        int nb_qbs;
};

class BinaryMatrixComputer : public MatrixComputer {
    public:
        BinaryMatrixComputer();
        BinaryMatrixComputer(int nb_gates, int nb_qbs);
        std::vector<std::vector<Eigen::MatrixXcd>> tree;
        void updateMatrix(int position, std::vector<std::shared_ptr<Gate>> list_gates);
        void calculateMatrix(std::vector<std::shared_ptr<Gate>> list_gates);
        void initializeTree(int nb_gates, int n_qubits);
        Eigen::MatrixXcd getMatrix();
        int nb_gates;
        int nb_qbs;
};
    

#endif