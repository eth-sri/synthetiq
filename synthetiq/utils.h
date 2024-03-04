#ifndef DEF_UTILS
#define DEF_UTILS

#include <Eigen/Dense>
#include <vector>

typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> BoolMatrix;

class Utils {
    public:
        static std::vector<int> dec2bin(int number, int size);
        static int bin2dec(std::vector<int> bin);
        static std::vector<int> getRowChangeQubitIndependent(std::vector<int> qbs, int nb_qbs);
        static Eigen::MatrixXcd increaseQubits(Eigen::MatrixXcd& matrix, int nb_qbs);
        static Eigen::MatrixXcd changeQubits(Eigen::MatrixXcd& matrix, std::vector<int> qbs_change);
        static bool matricesEqual(Eigen::MatrixXcd& matrix1, Eigen::MatrixXcd& matrix2);
        static std::vector<std::string> readQasmLineBase(std::string line);
        static std::string readNameQasmLine(std::string line);
        static int readNbQbsQasmLine(std::string line);
        static std::vector<int> readActingQubitsQasmLine(std::string line);
        static std::complex<double> traceConjugateProduct(Eigen::MatrixXcd& matrix1, Eigen::MatrixXcd& matrix2);
    private:
        static constexpr double epsilon = 1e-6;
};
    
#endif