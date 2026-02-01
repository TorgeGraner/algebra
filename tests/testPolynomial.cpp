#include <iostream>
#include <string>

#include "util/helper.hpp"

#include "datastructure/integer.hpp"
#include "datastructure/modp.hpp"
#include "datastructure/fractionalize.hpp"
#include "datastructure/epsFloat.hpp"

#include "datastructure/polynomial.hpp"

const int eps = 10000000;

template<typename R> int testDivision();

int main() {
    int numErr, numErrTotal = 0;

    numErr = 0;
    numErr += testDivision<ModP<2>>();
    numErr += testDivision<ModP<5>>();
    numErr += testDivision<Integer>();
    numErr += testDivision<Fractionalize<Integer>>();
    numErr += testDivision<EpsFloat<eps>>();

    numErrTotal += numErr;
    std::cout << "Encountered " << numErr << " errors testing polynomial division.\n";

    std::cout << "Finished testPolynomial with " << numErrTotal << " errors.\n";
    return numErrTotal;
}

template<typename R>
int testDivision() {
    std::vector<std::string> lhsStrArr, rhsStrArr, gTruthStrArr;

    lhsStrArr.emplace_back("-360 -450 -413 -618 -583 -216 -126 -163 42 7");
    lhsStrArr.emplace_back("-360 -450 -413 -618 -583 -216 -126 -163 42 7");
    lhsStrArr.emplace_back("-243 405 -270 90 -15 1");

    rhsStrArr.emplace_back("12 5 0 7");
    rhsStrArr.emplace_back("1 1 1 1 1");
    rhsStrArr.emplace_back("-3 1");

    gTruthStrArr.emplace_back("-30 -25 -24 -24 -24 6 1");
    gTruthStrArr.emplace_back("-360 -90 37 -205 35 7");
    gTruthStrArr.emplace_back("81 -108 54 -12 1");

    int cnt = -1;
    int numErr = 0;
    for (auto x : lhsStrArr) {
        ++cnt;
        // Initialization
        Polynomial<R> lhs, rhs, gTruth;
        Polynomial<R>::parseToPolynomial(lhsStrArr[cnt], lhs);
        Polynomial<R>::parseToPolynomial(rhsStrArr[cnt], rhs);
        Polynomial<R>::parseToPolynomial(gTruthStrArr[cnt], gTruth);

        try {
            // Function call
            Polynomial<R> res1 = lhs / rhs;

            lhs /= rhs;
            Polynomial<R> res2 = lhs;

            // Assertion
            numErr += util::assertEq(res1, gTruth, lhs, "polynomial division (out of place)");
            numErr += util::assertEq(res2, gTruth, lhs, "polynomial division (in place)");
        } catch (std::invalid_argument const& ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            ++numErr;
        }
    }
    if (cnt == -1) std::cout << "Warning: Did not test.\n";
    return numErr;
}
