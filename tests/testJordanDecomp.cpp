#include <iostream>
#include <string>

#include "util/helper.hpp"

#include "datastructure/integer.hpp"
#include "datastructure/modp.hpp"
#include "datastructure/fractionalize.hpp"
#include "datastructure/epsFloat.hpp"

#include "datastructure/matrix.hpp"

#include "algorithm/jordanDecomposition.hpp"

const int eps = 10000000;

template<typename R> int testJordanDecomposition();

int main() {
    int numErr, numErrTotal = 0;

    numErr = 0;
    numErr += testJordanDecomposition<Integer>();
    numErr += testJordanDecomposition<Fractionalize<Integer>>();
    numErr += testJordanDecomposition<ModP<5>>();
    numErr += testJordanDecomposition<ModP<2>>();
    numErr += testJordanDecomposition<EpsFloat<eps>>();

    numErrTotal += numErr;
    std::cout << "Encountered " << numErr << " errors testing jordan decomposition.\n";
    
    std::cout << "Finished testing decomposition with " << numErrTotal << " errors.\n";
    return numErrTotal;
}

template<typename R>
int testJordanDecomposition() {
    std::vector<std::string> lhsStrArr;

    lhsStrArr.emplace_back("-4 -4 8 -4 2 4 8 4 -4");
    lhsStrArr.emplace_back("1 0 1 0 1 1 1 1 0");
    lhsStrArr.emplace_back("1 2 0 0 3 0 2 -4 2");
    lhsStrArr.emplace_back("1 3 7 2 2 7 2 3 6");
    lhsStrArr.emplace_back("1 1 0 1 1 0 0 0 1");
    lhsStrArr.emplace_back("4 -9 6 12 0 -1 4 6 2 -11 8 16 -1 3 0 -1");
    lhsStrArr.emplace_back("4 0 0 0 0 4 0 0 1 -4 -4 0 -1 2 0 -4");
    lhsStrArr.emplace_back("0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0");
    lhsStrArr.emplace_back("1 -1 1 -1 0 -2 0 -3 2 3 0 3 1 5 -1 6");
    lhsStrArr.emplace_back("5 4 2 1 0 1 -1 -1 -1 -1 3 0 1 1 -1 2");
    lhsStrArr.emplace_back("1 2 2 1 -4 -1 0 1 -2 2 -3 2 6 2 -5 -3 -1 3 2 -2 -4 -1 4 0 0");
    lhsStrArr.emplace_back("1 2 3 4 5 0 1 2 3 4 0 0 1 2 3 0 0 0 1 2 0 0 0 0 1");
    lhsStrArr.emplace_back("2 1 0 0 0 0 2 0 0 0 0 0 2 0 0 0 0 0 2 1 0 0 0 0 2");
    lhsStrArr.emplace_back("4 -4 9 7 11 1 0 4 4 6 0 0 1 -1 -1 0 0 1 2 0 0 0 0 1 3");
    lhsStrArr.emplace_back("25 -16 30 -44 -12 13 -7 18 -26 -6 -18 12 -21 36 12 -9 6 -12 21 6 11 -8 15 -22 -3");
    lhsStrArr.emplace_back("2 1 0 0 0 0 -1 1 1 0 0 0 1 0 0 0 0 0 -1 0 2 3 1 0 1 0 -2 -1 1 0 -1 0 2 1 3 4");

    int cnt = -1;
    int numErr = 0;
    for (auto x : lhsStrArr) {
        ++cnt;
        // Initialization
        Matrix<R> lhs;
        Matrix<R>::parseToMatrix(lhsStrArr[cnt], lhs);

        try {
            // Function call
            Matrix<R> base = jordanDecomposition(lhs);
            Matrix<R> inv;
            R scale = base.inverse(inv);
            Matrix<R> decomp = inv * lhs * base;

            bool upperTriDiagonal = true;
            for (int i = 0; i < lhs.getN(); ++i) {
                for (int j = 0; j < lhs.getM(); ++j) {
                    int diff = i - j;
                    // diff < 0 above diagonal
                    if ((diff < -1 || diff > 0) && decomp(i, j) != 0) upperTriDiagonal = false;
                }
            }

            // Assertion
            numErr += util::assertEq(upperTriDiagonal, true, lhs, "jordan decomposition");
        } catch (std::invalid_argument const& ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            ++numErr;
        }
    }
    if (cnt == -1) std::cout << "Warning: Did not test.\n";
    return numErr;
}