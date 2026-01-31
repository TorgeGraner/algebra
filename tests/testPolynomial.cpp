#include <iostream>
#include <string>

#include "util/helper.hpp"

#include "datastructure/integer.hpp"
#include "datastructure/modp.hpp"
#include "datastructure/fractionalize.hpp"

#include "datastructure/polynomial.hpp"

template<typename R> int testDivision();

int main() {
    int numErrTotal = 0;
    int numErr = 0;

    numErr += testDivision<Integer>();
    numErr += testDivision<Fractionalize<Integer>>();

    std::cout << numErr << " errors occurred\n";
    numErrTotal += numErr;
    numErr = 0;
    return numErrTotal;
}

template<typename R>
int testDivision() {
    // Not that for mod p coefficients, the determinant has to be nonzero mod p
    std::string dividends[] = {
        "-360 -450 -413 -618 -583 -216 -126 -163 42 7",
        "-360 -450 -413 -618 -583 -216 -126 -163 42 7",
        "-243 405 -270 90 -15 1"
    };
    std::string divisors[] = {
        "12 5 0 7",
        "1 1 1 1 1",
        "-3 1 "
    };
    std::string results[] = {
        "-30 -25 -24 -24 -24 6 1",
        "-360 -90 37 -205 35 7",
        "81 -108 54 -12 1"
    };
    int cnt = 0;
    for (std::string str : dividends) {
        Polynomial<R> divisor, dividend, result;
        Polynomial<R>::parseToPolynomial(dividends[cnt], dividend);
        Polynomial<R>::parseToPolynomial(divisors[cnt], divisor);
        Polynomial<R>::parseToPolynomial(results[cnt], result);
        try {
            std::cout << dividend / divisor << " " << result << "\n";
            if (dividend / divisor != result) return EXIT_FAILURE;
        } catch (std::invalid_argument const& ex) {
            std::cout << ex.what() << "\n";
            return EXIT_FAILURE;
        }
        ++cnt;
    }
    return EXIT_SUCCESS;
}
