#pragma once
#include <vector>
#include <sstream>

#include "util/helper.hpp"

#include "datastructure/integer.hpp"
#include "datastructure/fractionalize.hpp"
#include "datastructure/matrix.hpp"
#include "datastructure/polynomial.hpp"

#include "algorithm/matrix_operations.hpp"
#include "algorithm/jordan_decomposition.hpp"

template<typename R> Matrix<R> parseToMatrix(std::string);
Matrix<int> fullBipartite(int n, int m);

int main() {
    while(true) {
        std::string line;
        std::string str = "";
        std::cout << "Input Matrix, enter with line break, exit with empty matrix:\n"; 
        do {
            getline(std::cin, line);
            str += line + " ";
        } while (!line.empty());

        auto fMat = parseToMatrix<Fractionalize<Integer>>(str);
        auto iMat = parseToMatrix<Integer>(str);

        if (fMat.getM() == 0) break;

        std::cout << "Decomposing matrix: " << std::endl << fMat << std::endl;

        std::cout << "Decomposing with integer coefficients:" << std::endl;
        auto iBase = decompose(iMat);
        auto iInv = matOps::inverse(iBase);
        std::cout << "Base:" << std::endl << iBase << std::endl << "Inverse: " << std::endl << iInv << std::endl << "Decomposition:" << std::endl << iInv * iMat * iBase << std::endl;

        std::cout << "Decomposing with rational coefficients:" << std::endl;
        auto fBase = decompose(fMat);
        auto fInv = matOps::inverse(fBase);
        std::cout << "Base:" << std::endl << fBase << std::endl << "Inverse: " << std::endl << fInv << std::endl << "Decomposition:" << std::endl << fInv * fMat * fBase << std::endl;
    }
    std::cout << "Finished" << std::endl;

    util::allocate<Integer>(0, true);
    util::allocate<Fractionalize<Integer>>(0, true);

    util::deallocate<Integer>(nullptr, true);
    util::deallocate<Fractionalize<Integer>>(nullptr, true);

    std::string line;
    getline(std::cin, line);
    return 1;
}

template<typename R>
Matrix<R> parseToMatrix(std::string str) {
    // Get dimension
    int dim_sq = 0;
    int num = 0;
    int cnt = 0;
    std::istringstream is(str);
    while (is >> num) dim_sq++;

    // Get values
    R* values = util::allocate<R>(dim_sq);
    is = std::istringstream (str);
    while (is >> num) values[cnt++] = num;

    int dim = (int) std::sqrt(dim_sq);
    Matrix<R> result(values, dim, dim);
    util::deallocate(values);
    return result;
}

Matrix<int> fullBipartite(int n, int m) {
    int k = n + m;
    int* values = util::allocate<int>(k * k);
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            values[k * i + j] = (((j < n && i < n) || (j >= n && i >= n)) ? 0 : 1);
        }
    }
    Matrix<int> ret(values, k, k);
    util::deallocate(values);
    return ret;
}