#include <iostream>
#include <string>

#include "util/helper.hpp"

#include "datastructure/integer.hpp"
#include "datastructure/modp.hpp"
#include "datastructure/fractionalize.hpp"

#include "datastructure/matrix.hpp"

#define deb = false;

template<typename R> int testInverse();

template<typename R> int testMult();

int main() {
    int numErrTotal = 0;

    numErrTotal += testMult<Integer>();

    std::cout << "Finished with " << numErrTotal << " errors.\n";

    return numErrTotal;
}

template<typename R>
int testMult() {
    std::string lhsStrArr[] = {
        "-4 -4 8 -4 2 4 8 4 -4",
    };
    std::string rhsStrArr[] = {
        "-4 -4 8 -4 2 4 8 4 -4",
    };
    std::string gTruthStrArr[] = {
        "96 40 -80 40 36 -40 -80 -40 96",
    };
    int dims[][3] = {
        {3, 3, 3},
    };
    int cnt = -1;
    int numErr = 0;
    for (auto x : lhsStrArr) {
        ++cnt;
        Matrix<R> lhs, rhs, gTruth;
        if (!Matrix<R>::parseToMatrix(lhsStrArr[cnt], lhs, dims[cnt][0]) || !Matrix<R>::parseToMatrix(rhsStrArr[cnt], rhs, dims[cnt][1]) || !Matrix<R>::parseToMatrix(gTruthStrArr[cnt], gTruth, dims[cnt][2])) {
            ++numErr;
            continue;
        }
        lhs *= rhs;
        numErr += (lhs != gTruth);
    }
    return numErr;
}


template<typename R>
int testInverse() {
    std::string sArr[] = {
        "1 -1 1 -1 0 -2 0 -3 2 3 0 3 1 5 -1 6"
    };
    
    for (std::string matStr : sArr) {
        try {
            Matrix<R> mat;
            if(!Matrix<R>::parseToMatrix(matStr, mat)) return EXIT_FAILURE;
            // Function call
            Matrix<R> inv;
            R scale = mat.inverse(inv);
            if (mat * inv != scale) return EXIT_FAILURE;
            if (inv * mat != scale) return EXIT_FAILURE;
            std::cout << "No error." << std::endl;
        } catch (std::invalid_argument ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
