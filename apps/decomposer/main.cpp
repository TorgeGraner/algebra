#pragma once
#include <vector>
#include <sstream>
#include <chrono>
#include <cstdlib>

#include "util/helper.hpp"

#include "datastructure/integer.hpp"
#include "datastructure/fractionalize.hpp"
#include "datastructure/matrix.hpp"
#include "datastructure/polynomial.hpp"
#include "datastructure/modp.hpp"

#include "algorithm/matrix_operations.hpp"
#include "algorithm/jordan_decomposition.hpp"

const int p = 2;

template<typename R> bool parseToMatrix(std::string, Matrix<R>&);
template<typename R> void decomp(Matrix<R>&);
Matrix<int> fullBipartite(int n, int m);

void account() {
    #ifndef NDEBUG
        std::cout << "Accounting:" << std::endl;

        util::allocate<Integer>(0, true);
        util::allocate<Fractionalize<Integer>>(0, true);
        util::allocate<ModP<p>>(0, true);

        util::deallocate<Integer>(nullptr, true);
        util::deallocate<Fractionalize<Integer>>(nullptr, true);
        util::deallocate<ModP<p>>(0, true);
    #endif
}

int main() {
    if (std::atexit(account)) return EXIT_FAILURE;

    while(true) {
        std::string line;
        std::string str = "";
        std::cout << "Input Matrix, enter with line break, exit with empty matrix:\n"; 
        do {
            getline(std::cin, line);
            str += line + " ";
        } while (!line.empty());

        Matrix<Integer> iMat;
        Matrix<Fractionalize<Integer>> fMat;
        Matrix<ModP<p>> mMat;

        if (!parseToMatrix<Fractionalize<Integer>>(str, fMat) || !parseToMatrix<Integer>(str, iMat) || !parseToMatrix<ModP<p>>(str, mMat)) {
            std::cout << "Erroneous matrix entered"<< std::endl;
            continue;
        }

        if (fMat.getM() == 0) return EXIT_SUCCESS;

        decomp(iMat);
        decomp(fMat);
        decomp(mMat);
    }
    
}

template<typename R>
void decomp(Matrix<R>& mat) {
    std::cout << "Decomposing with " << typeid(R).name() << " coefficients:" << std::endl;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    auto base = decompose(mat);
    auto inv = matOps::inverse(base);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Base:" << std::endl << base << "Inverse: " << std::endl << inv << "Decomposition:" << std::endl << inv * mat * base;
    #ifndef NDEBUG
        std::cout << "Required " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " microseconds" << std::endl;
    #endif
    std::cout << std::endl;
}

template<typename R>
bool parseToMatrix(std::string str, Matrix<R>& result) {
    // Get dimension
    int dim_sq = 0;
    int num = 0;
    int cnt = 0;
    std::istringstream is(str);
    while (is >> num) dim_sq++;
    int dim = (int) std::sqrt(dim_sq);
    if (dim * dim != dim_sq) return false;

    // Get values
    R* values = util::allocate<R>(dim_sq);
    is = std::istringstream (str);
    while (is >> num) values[cnt++] = num;

    result = Matrix<R>(values, dim, dim);
    util::deallocate(values);
    return true;
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