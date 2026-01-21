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
#include "datastructure/root.hpp"

#include "algorithm/matrix_operations.hpp"
#include "algorithm/jordan_decomposition.hpp"

// Must be prime
const int p = 2;

template<typename R> bool parseToMatrix(std::string, Matrix<R>&);
template<typename R> void decomp(Matrix<R>&);
template<typename R> Matrix<R> fullBipartite(int n, int m);

void account() { 
    std::cout << "There were " << util::numAlloc - util::numDealloc << " unaccounted arrays.\n"; 
}

int main() {
    #ifndef NDEBUG
        if (std::atexit(account)) return EXIT_FAILURE;
    #endif

    while(true) {
        std::string line;
        std::string str = "";
        std::cout << "Input Matrix or exit by typing exit:\n"; 
        do {
            getline(std::cin, line);
            str += line + " ";
        } while (!line.empty());

        Matrix<Integer> iMat;
        Matrix<Fractionalize<Integer>> fMat;
        Matrix<ModP<p>> mMat;

        if (str.find("exit") != std::string::npos) {
            return EXIT_SUCCESS;
        }
        if (str.find_first_not_of(" ") == std::string::npos) {
            std::cout << "Please enter a matrix\n";
            continue;
        }
        if (!parseToMatrix<Fractionalize<Integer>>(str, fMat) || !parseToMatrix<Integer>(str, iMat) || !parseToMatrix<ModP<p>>(str, mMat)) {
            std::cout << "Erroneous matrix entered\n";
            continue;
        }

        decomp(iMat);
        decomp(fMat);
        decomp(mMat);
    }
}

template<typename R>
void decomp(Matrix<R>& mat) {
    std::cout << "Decomposing with " << typeid(R).name() << " coefficients:\n";

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    Matrix<R> base = decompose(mat);
    Matrix<R> inv;
    R scale = matOps::inverse(base, inv);
    Matrix<R> decomposition = inv * mat * base;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Base:\n" << base;
    if (scale != 1) {
        std::cout << "Inverse scaled by " << scale << ":\n" << inv << "Decomposition scaled by " << scale << ":\n" << decomposition;
    } else {
        std::cout << "Inverse:\n" << inv << "Decomposition:\n" << decomposition;
    }

    #ifndef NDEBUG
        std::cout << "Required " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " microseconds\n";
    #endif
    std::cout << "\n";
}

/*
* @brief Returns a Matrix with coefficients in R from a string of integers. The number of integers in the string must be a perfect square.
*/
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

/*
* @brief Returns the adjacency matrix of the full bipartite graph
*/
template<typename R>
Matrix<R> fullBipartite(int n, int m) {
    int k = n + m;
    R* values = util::allocate<R>(k * k);
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            values[k * i + j] = (((j < n && i < n) || (j >= n && i >= n)) ? 0 : 1);
        }
    }
    Matrix<R> ret(values, k, k);
    util::deallocate(values);
    return ret;
}