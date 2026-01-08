#pragma once
#include <vector>

#include "util/helper.hpp"

#include "datastructure/integer.hpp"
#include "datastructure/fractionalize.hpp"
#include "datastructure/matrix.hpp"
#include "datastructure/polynomial.hpp"

#include "algorithm/matrix_operations.hpp"
#include "algorithm/jordan_decomposition.hpp"


Matrix<int> fullBipartite(int n, int m);

int main() {
    int values[] = {

        // Diagonalizable


        //-4, -4, 8, -4, 2, 4, 8, 4, -4
        //1, 0, 1, 0, 1, 1, 1, 1, 0
        //1, 2, 0, 0, 3, 0, 2, -4, 2
        //1, 3, 7, 2, 2, 7, 2, 3, 6

        //4, 0, 0, 0, 0, 4, 0, 0, 1, -4, -4, 0, -1, 2, 0, -4
        //4, -9, 6, 12, 0, -1, 4, 6, 2, -11, 8, 16, -1, 3, 0, -1

        //0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0

        // Trigonalizable:

         //1, -1, 1, -1, 0, -2, 0, -3, 2, 3, 0, 3, 1, 5, -1, 6
         //5, 4, 2, 1, 0, 1, -1, -1, -1, -1, 3, 0, 1, 1, -1, 2

         //1, 2, 2, 1, -4, -1, 0, 1, -2, 2, -3, 2, 6, 2, -5, -3, -1, 3, 2, -2, -4, -1, 4, 0, 0
         //1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 0, 0, 1, 2, 3, 0, 0, 0, 1, 2, 0, 0, 0, 0, 1
         //2, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 2
         //4, -4, 9, 7, 11, 1, 0, 4, 4, 6, 0, 0, 1, -1, -1, 0, 0, 1, 2, 0, 0, 0, 0, 1, 3

         25, -16, 30, -44, -12, 13, -7, 18, -26, -6, -18, 12, -21, 36, 12, -9, 6, -12, 21, 6, 11, -8, 15, -22, -3

         //2, 1, 0, 0, 0, 0, -1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 2, 3, 1, 0, 1, 0, -2, -1, 1, 0, -1, 0, 2, 1, 3, 4

         // no rational roots

         //2, 1, 1, 1, -2, 1, 1, 1, 0, -1, 1, 0, 2, 0, -1, 1, 0, 1, 2, -2, 1, 0, 1, 0, 0
         //2,	1,	0,  0,	0,	0,	0,	0, 1,	2,	1,	0,	0,	0,	0,	0, 0,	1,	2,	1,	0,	0,	0,	0, 0,	0,	1,	2,	1,	0,	0,	0, 0,	0,	0,	1,	2,	1,	0,	0, 0,	0,	0,	0,	1,	2,	1,	0, 0,	0,	0,	0,	0,	1,	2,	1, 0,	0,	0,	0,	0,	0,	1,	2

         // Test spanBasis

         //3, 0, 0, 0, 1, 0, 0, 0, 5, 0, 2, 0, 0, 0, 0, 7, 4, 0
         //1, 2, 2, 7, 1, 5, 6, 2, 3, 3, 4, 5, 3, 3, 6, 7

    };

    int l = sizeof(values) / sizeof(int);
    int s = (int) sqrt(l);

    Fractionalize<Integer>* fValues = util::allocate<Fractionalize<Integer>>(l);
    Integer* iValues = util::allocate<Integer>(l);

    for (int i = 0; i < l; ++i) {
        fValues[i] = Integer(values[i]);
        iValues[i] = Integer(values[i]);
    }

    Matrix<Fractionalize<Integer>> fMat(fValues, s, s);
    Matrix<Integer> iMat(iValues, s, s);
    util::deallocate(fValues);
    util::deallocate(iValues);

    std::cout << "To decompose: " << std::endl << fMat << std::endl;
    try {
        std::cout << "Decomposing with integer coefficients:" << std::endl;
        auto iBase = decompose(iMat);
        auto iInv = matOps::inverse(iBase);
        std::cout << "Base:" << std::endl << iBase << std::endl << "Inverse: " << std::endl << iInv << std::endl << "Decomposition:" << std::endl << iInv * iMat * iBase << std::endl;

        std::cout << "Decomposing with rational coefficients:" << std::endl;
        auto fBase = decompose(fMat);
        auto fInv = matOps::inverse(fBase);
        std::cout << "Base:" << std::endl << fBase << std::endl << "Inverse: " << std::endl << fInv << std::endl << "Decomposition:" << std::endl << fInv * fMat * fBase << std::endl;

    }
    catch (std::invalid_argument const& ex) {
        std::cout << "Invalid argument exception: " << ex.what() << std::endl;
        return -1;
    }
    catch (std::out_of_range const& ex) {
        std::cout << "Out of range exception: " << ex.what() << std::endl;
        return -1;
    }

    util::allocate<Integer>(0, true);
    util::allocate<Fractionalize<Integer>>(0, true);

    util::deallocate<Integer>(nullptr, true);
    util::deallocate<Fractionalize<Integer>>(nullptr, true);
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
