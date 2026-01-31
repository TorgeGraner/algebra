#include <iostream>
#include <string>

#include "util/helper.hpp"

#include "datastructure/integer.hpp"
#include "datastructure/modp.hpp"
#include "datastructure/fractionalize.hpp"

#include "datastructure/matrix.hpp"

#include "algorithm/jordan_decomposition.hpp"

template<typename R> int testJordanDecomp();

int main() {
    int numErrTotal = 0;
    int numErr = 0;

    // Test Inverse
    numErr += testJordanDecomp<Integer>();
    numErr += testJordanDecomp<Fractionalize<Integer>>();
    numErr += testJordanDecomp<ModP<5>>();
    numErr += testJordanDecomp<ModP<2>>();

    numErrTotal += numErr;
    std::cout << "Finished with " << numErrTotal << " errors.\n";
    return numErrTotal;
}

template<typename R>
int testJordanDecomp() {
    int numErr = 0;
    std::string sArr[] = {
        "-4 -4 8 -4 2 4 8 4 -4",
        "1 0 1 0 1 1 1 1 0",
        "1 2 0 0 3 0 2 -4 2",
        "1 3 7 2 2 7 2 3 6",
        "1 1 0 1 1 0 0 0 1",
        "4 -9 6 12 0 -1 4 6 2 -11 8 16 -1 3 0 -1",
        "4 0 0 0 0 4 0 0 1 -4 -4 0 -1 2 0 -4",
        "0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 ",
        "1 -1 1 -1 0 -2 0 -3 2 3 0 3 1 5 -1 6",
        "5 4 2 1 0 1 -1 -1 -1 -1 3 0 1 1 -1 2",
        "1 2 2 1 -4 -1 0 1 -2 2 -3 2 6 2 -5 -3 -1 3 2 -2 -4 -1 4 0 0",
        "1 2 3 4 5 0 1 2 3 4 0 0 1 2 3 0 0 0 1 2 0 0 0 0 1",
        "2 1 0 0 0 0 2 0 0 0 0 0 2 0 0 0 0 0 2 1 0 0 0 0 2",
        "4 -4 9 7 11 1 0 4 4 6 0 0 1 -1 -1 0 0 1 2 0 0 0 0 1 3",
        "25 -16 30 -44 -12 13 -7 18 -26 -6 -18 12 -21 36 12 -9 6 -12 21 6 11 -8 15 -22 -3",
        "2 1 0 0 0 0 -1 1 1 0 0 0 1 0 0 0 0 0 -1 0 2 3 1 0 1 0 -2 -1 1 0 -1 0 2 1 3 4",
    };
    
    for (std::string matStr : sArr) {
        Matrix<R> mat, base, inv, decomp;
        R scale;
        try {
            Matrix<R> mat;
            if(!Matrix<R>::parseToMatrix(matStr, mat)) return EXIT_FAILURE;
            // Function call
            base = decompose(mat);
            scale = base.inverse(inv);
            decomp = inv * mat * base;
            // Condition
            // Assert the matrix is in jordan normal form
            for (int i = 0; i < mat.getN(); ++i) {
                for (int j = 0; j < mat.getM(); ++j) {
                    int diff = i - j;
                    if ((diff > 0 || diff < -1) && decomp(i, j) != 0) {
                        std::cerr << "Error using matrix:\n" << mat << "\n";
                        std::cerr << "Base:\n" << base << "Inverse:\n" << inv << "Decomposition:\n" << decomp << "\n";
                        i = mat.getN();
                        break;
                        ++numErr;
                    }
                }
            }
        } catch (std::invalid_argument ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            std::cout << "Base:\n" << base << "Inverse:\n" << inv << "Decomposition:\n" ;
            ++numErr;
        }
    }
    return numErr;
}

/*
// Diagonalizable (over Q, possibly not diagonalizable over F_p, especially F_2)
-4 -4 8 -4 2 4 8 4 -4
1 0 1 0 1 1 1 1 0
1 2 0 0 3 0 2 -4 2
1 3 7 2 2 7 2 3 6
1 1 0 1 1 0 0 0 1
4 -9 6 12 0 -1 4 6 2 -11 8 16 -1 3 0 -1
4 0 0 0 0 4 0 0 1 -4 -4 0 -1 2 0 -4
0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 

// Trigonalizable:
1 -1 1 -1 0 -2 0 -3 2 3 0 3 1 5 -1 6
5 4 2 1 0 1 -1 -1 -1 -1 3 0 1 1 -1 2
1 2 2 1 -4 -1 0 1 -2 2 -3 2 6 2 -5 -3 -1 3 2 -2 -4 -1 4 0 0
1 2 3 4 5 0 1 2 3 4 0 0 1 2 3 0 0 0 1 2 0 0 0 0 1
2 1 0 0 0 0 2 0 0 0 0 0 2 0 0 0 0 0 2 1 0 0 0 0 2
4 -4 9 7 11 1 0 4 4 6 0 0 1 -1 -1 0 0 1 2 0 0 0 0 1 3
25 -16 30 -44 -12 13 -7 18 -26 -6 -18 12 -21 36 12 -9 6 -12 21 6 11 -8 15 -22 -3
2 1 0 0 0 0 -1 1 1 0 0 0 1 0 0 0 0 0 -1 0 2 3 1 0 1 0 -2 -1 1 0 -1 0 2 1 3 4

// Not enough integer roots for a full di-/triagonalization (using this algorithm)
2 1 1 1 -2 1 1 1 0 -1 1 0 2 0 -1 1 0 1 2 -2 1 0 1 0 0
2 1 0 0 0 0 0 0 1 2 1 0 0 0 0 0 0 1 2 1 0 0 0 0 0 0 1 2 1 0 0 0 0 0 0 1 2 1 0 0 0 0 0 0 1 2 1 0 0 0 0 0 0 1 2 1 0 0 0 0 0 0 1 2

*/