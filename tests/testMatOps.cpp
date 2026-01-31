#include <iostream>
#include <string>

#include "util/helper.hpp"

#include "datastructure/integer.hpp"
#include "datastructure/modp.hpp"
#include "datastructure/fractionalize.hpp"

#include "datastructure/matrix.hpp"

#include "algorithm/matrix_operations.hpp"

template<typename R> int testInverse();
template<typename R> int testRank();
template<typename R> int testDeterminant();
template<typename R> int testKernelBasis();
template<typename R> int testCompleteBasis();

/*
int testCharPoly();
int testCharPolyNaive();
*/

int main() {
    int numErrTotal = 0;
    int numErr = 0;

    // Test Inverse
    numErr += testInverse<Integer>();
    numErr += testInverse<Fractionalize<Integer>>();
    numErr += testInverse<ModP<5>>();
    numErr += testInverse<ModP<2>>();
    if(numErr > 0) std::cout << "Finished testing inverse() with " << numErr << " errors.\n";
    numErrTotal += numErr;
    numErr = 0;

    numErr += testRank<Integer>();
    numErr += testRank<Fractionalize<Integer>>();
    numErr += testRank<ModP<5>>();
    numErr += testRank<ModP<2>>();

    if(numErr > 0) std::cout << "Finished testing rank() with " << numErr << " errors.\n";
    numErrTotal += numErr;
    numErr = 0;

    numErr += testDeterminant<Integer>();
    numErr += testDeterminant<Fractionalize<Integer>>();
    numErr += testDeterminant<ModP<5>>();
    numErr += testDeterminant<ModP<2>>();

    if(numErr > 0) std::cout << "Finished testing determinant() with " << numErr << " errors.\n";
    numErrTotal += numErr;
    numErr = 0;

    numErr += testKernelBasis<Integer>();
    numErr += testKernelBasis<Fractionalize<Integer>>();
    numErr += testKernelBasis<ModP<5>>();
    numErr += testKernelBasis<ModP<2>>();

    if(numErr > 0) std::cout << "Finished testing kernelBasis() with " << numErr << " errors.\n";
    numErrTotal += numErr;
    numErr = 0;
    return numErrTotal;
}

template<typename R>
int testInverse() {
    // Not that for mod p coefficients, the determinant has to be nonzero mod p
    std::string sArr[] = {
        "1 4 0 2 3 6 0 6 7", // det -71
        "0 0 1 0 1 0 1 0 0", // det -1
        "1 0 2 0 1 0 1 0 3", // det 1
    };

    for (std::string matStr : sArr) {
        Matrix<R> mat;
        try {
            if(!Matrix<R>::parseToMatrix(matStr, mat)) return EXIT_FAILURE;
            Matrix<R> inv;
            // Function call
            R scale = matOps::inverse(mat, inv);
            // Condition
            if(mat * inv != scale || inv * mat != scale) return EXIT_FAILURE;
        } catch (std::invalid_argument ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

template<typename R>
int testRank() {
    std::string sArr[] = {
        "3 1 1 2 1",
        "1 2 3 2 1",
        "0 0 0 0",
        "0 3 1 0 3 3 3 3 1 3 3 0 0 1 4 4",
        "0 3 0 0 1 1 3 3 1 1 1 0 0 1 0 0"
    };
    int dims[] = {
        5,
        1,
        2,
        4,
        4
    };
    int ranks[] = {
        1,
        1,
        0,
        4,
        3
    };

    int cnt = 0;
    for (std::string matStr : sArr) {
        Matrix<R> mat;
        try {
            if(!Matrix<R>::parseToMatrix(matStr, mat, dims[cnt])) return EXIT_FAILURE;
            // Function call
            int rk = matOps::rank(mat);
            // Condition
            if(rk != ranks[cnt]) {
                std::cerr << "Error using " << typeid(R).name() << " coefficients. Rank() returned " << rk << " when it was " << ranks[cnt] << "\n";
                return EXIT_FAILURE;
            }
        } catch (std::invalid_argument ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            return EXIT_FAILURE;
        }
        ++cnt;
    }
    return EXIT_SUCCESS;
}

template<typename R>
int testDeterminant() {
    // Not that for mod p coefficients, the determinant has to be nonzero mod p
    std::string sArr[] = {
        "1 4 0 2 3 6 0 6 7", // det -71
        "0 0 1 0 1 0 1 0 0", // det -1
        "1 0 2 0 1 0 1 0 3", // det 1
        "0"
    };
    int dets [] = {
        -71,
        -1,
        1,
        0
    };

    int cnt = 0;
    for (std::string matStr : sArr) {
        Matrix<R> mat;
        try {
            if(!Matrix<R>::parseToMatrix(matStr, mat)) return EXIT_FAILURE;
            // Function call
            R det = matOps::determinant(mat);
            // Condition
            if(det != dets[cnt]) return EXIT_FAILURE;
        } catch (std::invalid_argument ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            return EXIT_FAILURE;
        }
        ++cnt;
    }
    return EXIT_SUCCESS;
}

template<typename R>
int testKernelBasis() {
    std::string sArr[] = {
        "1 4 0 2 3 6 0 6 7", // det -71
        "0 0 1 0 1 0 1 0 0", // det -1
        "1 0 2 0 1 0 1 0 3", // det 1
        "0"
    };
    int cnt = 0;
    for (std::string matStr : sArr) {
        Matrix<R> mat;
        try {
            if(!Matrix<R>::parseToMatrix(matStr, mat)) return EXIT_FAILURE;
            // Function call
            Matrix<R> kBase= matOps::kernelBasis(mat);
            // Condition
            if (mat * kBase != 0 && matOps::rank(mat) == mat.getN() - matOps::rank(kBase)) return EXIT_FAILURE;
        } catch (std::invalid_argument ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            return EXIT_FAILURE;
        }
        ++cnt;
    }
    return EXIT_SUCCESS;
}

template<typename R>
int testCompleteBasis() {
    return EXIT_SUCCESS;
}

/*
testCharPoly(const Matrix<R>&);
testCharPolyNaive(const Matrix<R>&);
*/
