#include <iostream>
#include <string>

#include "util/helper.hpp"

#include "datastructure/integer.hpp"
#include "datastructure/modp.hpp"
#include "datastructure/epsFloat.hpp"
#include "datastructure/fractionalize.hpp"

#include "datastructure/matrix.hpp"
#include "datastructure/polynomial.hpp"

const int eps = 10000000;

template<typename R> int testMultiply();
template<typename R> int testRank();
template<typename R> int testDeterminant();
template<typename R> int testInverse();
template<typename R> int testNullspace();
template<typename R> int testComplete();
template<typename R> int testCharPoly();

int main() {

    int numErr, numErrTotal = 0;

    numErr = 0;
    numErr += testMultiply<Integer>();
    numErr += testMultiply<Fractionalize<Integer>>();
    numErr += testMultiply<ModP<3>>();
    numErr += testMultiply<ModP<5>>();
    numErr += testMultiply<EpsFloat<eps>>();
    
    numErrTotal += numErr;
    std::cout << "Encountered " << numErr << " errors testing matrix multiplication.\n";
    
    numErr = 0; 
    numErr += testInverse<Integer>();
    numErr += testInverse<Fractionalize<Integer>>();
    numErr += testInverse<ModP<3>>();
    numErr += testInverse<ModP<5>>();
    numErr += testInverse<EpsFloat<eps>>();

    numErrTotal += numErr;
    std::cout << "Encountered " << numErr << " errors testing matrix inversion.\n";

    numErr = 0; 
    numErr += testRank<Integer>();
    numErr += testRank<Fractionalize<Integer>>();
    numErr += testRank<ModP<2>>();
    numErr += testRank<ModP<5>>();
    numErr += testRank<EpsFloat<eps>>();

    numErrTotal += numErr;
    std::cout << "Encountered " << numErr << " errors testing matrix rank.\n";

    numErr = 0; 
    numErr += testDeterminant<Integer>();
    numErr += testDeterminant<Fractionalize<Integer>>();
    numErr += testDeterminant<ModP<2>>();
    numErr += testDeterminant<ModP<5>>();
    numErr += testDeterminant<EpsFloat<eps>>();

    numErrTotal += numErr;
    std::cout << "Encountered " << numErr << " errors testing matrix determinant.\n";

    numErr = 0; 
    numErr += testNullspace<Integer>();
    numErr += testNullspace<Fractionalize<Integer>>();
    numErr += testNullspace<ModP<7>>();
    numErr += testNullspace<ModP<5>>();
    numErr += testNullspace<EpsFloat<eps>>();

    numErrTotal += numErr;
    std::cout << "Encountered " << numErr << " errors testing nullspace.\n";
    
    numErr = 0; 
    numErr += testCharPoly<Integer>();
    numErr += testCharPoly<Fractionalize<Integer>>();
    numErr += testCharPoly<ModP<2>>();
    numErr += testCharPoly<ModP<5>>();
    numErr += testCharPoly<EpsFloat<eps>>();


    numErrTotal += numErr;
    std::cout << "Encountered " << numErr << " errors testing charpoly.\n";
    
    numErr = 0; 
    numErr += testComplete<Integer>();
    numErr += testComplete<Fractionalize<Integer>>();
    numErr += testComplete<ModP<2>>();
    numErr += testComplete<ModP<5>>();
    numErr += testComplete<EpsFloat<eps>>();

    numErrTotal += numErr;
    std::cout << "Encountered " << numErr << " errors testing complete.\n";
    std::cout << "Finished testing class Matrix with " << numErrTotal << " errors.\n";

    return numErrTotal;
}

template<typename R>
int testMultiply() {
    std::vector<std::string> lhsStrArr, rhsStrArr, gTruthStrArr;
    std::vector<int> lhsDimArr, rhsDimArr, gTruthDimArr;

    lhsStrArr.emplace_back("-4 -4 8 -4 2 4 8 4 -4");

    rhsStrArr.emplace_back("-4 -4 8 -4 2 4 8 4 -4");

    gTruthStrArr.emplace_back("96 40 -80 40 36 -40 -80 -40 96");
    
    lhsDimArr.emplace_back(3);

    rhsDimArr.emplace_back(3);
    
    gTruthDimArr.emplace_back(3);


    int cnt = -1;
    int numErr = 0;
    for (auto x : lhsStrArr) {
        ++cnt;
        // Initialization
        Matrix<R> lhs, rhs, gTruth;
        Matrix<R>::parseToMatrix(lhsStrArr[cnt],    lhs,    lhsDimArr[cnt]);
        Matrix<R>::parseToMatrix(rhsStrArr[cnt],    rhs,    rhsDimArr[cnt]);
        Matrix<R>::parseToMatrix(gTruthStrArr[cnt], gTruth, gTruthDimArr[cnt]);

        try {
            // Function call
            Matrix<R> res1 = lhs * rhs;
            lhs *= rhs;
            Matrix<R> res2 = lhs;
            // Assertion
            numErr += util::assertEq(res1, gTruth, lhs, "matrix multiplication (out of place)");
            numErr += util::assertEq(res2, gTruth, lhs, "matrix multiplication (in place)");
        } catch (std::invalid_argument const& ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            ++numErr;
        }
    }
    if (cnt == -1) std::cout << "Warning: Did not test.\n";
    return numErr;
}

template<typename R>
int testInverse() {
    std::vector<std::string> lhsStrArr;
    lhsStrArr.emplace_back("1 -1 1 -1 0 -2 0 -3 2 3 0 3 1 5 -1 6"); // det 2
    lhsStrArr.emplace_back("16 36 25 5 0 0 8 12 4 0 0 0 4 3 0 0 0 0 2 0 0 0 0 0 1");

    int cnt = -1;
    int numErr = 0;
    for (std::string str : lhsStrArr) {
        ++cnt;

        // Initialization
        Matrix<R> lhs;
        Matrix<R>::parseToMatrix(lhsStrArr[cnt], lhs);

        try {
            // Function call
            Matrix<R> inv;
            R scale = lhs.inverse(inv);
            Matrix<R> leftMult = lhs * inv;
            Matrix<R> rightMult = inv * lhs;
            // Assertion
            numErr += util::assertEq(leftMult, scale, lhs, "matrix inversion (left sided multiplication)");
            numErr += util::assertEq(rightMult, scale, lhs, "matrix inversion (right sided multiplication)");
        } catch (std::invalid_argument ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            ++numErr;
        }
    }
    if (cnt == -1) std::cout << "Warning: Did not test.\n";
    return numErr;
}

template<typename R>
int testRank() {
    std::vector<std::string> lhsStrArr;
    std::vector<int> lhsDimArr;
    std::vector<int> gTruthArr;

    lhsStrArr.emplace_back("3 1 1 2 1");
    lhsStrArr.emplace_back("1 2 3 2 1");
    lhsStrArr.emplace_back("0 0 0 0");
    lhsStrArr.emplace_back("0 3 1 0 3 3 3 3 1 3 3 0 0 1 4 4");
    lhsStrArr.emplace_back("0 3 0 0 1 1 3 3 1 1 1 0 0 1 0 0");

    lhsDimArr.emplace_back(5);
    lhsDimArr.emplace_back(1);
    lhsDimArr.emplace_back(2);
    lhsDimArr.emplace_back(4);
    lhsDimArr.emplace_back(4);

    gTruthArr.emplace_back(1);
    gTruthArr.emplace_back(1);
    gTruthArr.emplace_back(0);
    gTruthArr.emplace_back(4);
    gTruthArr.emplace_back(3);

    int cnt = -1;
    int numErr = 0;
    for (std::string str : lhsStrArr) {
        ++cnt; 
        // Initialization
        Matrix<R> lhs;
        Matrix<R>::parseToMatrix(lhsStrArr[cnt], lhs, lhsDimArr[cnt]);
        try {
            // Function call
            int rank = lhs.getRank();
            // Condition
            numErr += util::assertEq(rank, gTruthArr[cnt], lhs, "matrix rank");
        } catch (std::invalid_argument ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            ++numErr;
        }
    }
    if (cnt == -1) std::cout << "Warning: Did not test.\n";
    return numErr;
}


template<typename R>
int testDeterminant() {
    std::vector<std::string> lhsStrArr;
    std::vector<int> gTruthArr;

    lhsStrArr.emplace_back("1 4 0 2 3 6 0 6 7");
    lhsStrArr.emplace_back("0 0 1 0 1 0 1 0 0");
    lhsStrArr.emplace_back("1 0 2 0 1 0 1 0 3");
    lhsStrArr.emplace_back("0");

    gTruthArr.emplace_back(-71);
    gTruthArr.emplace_back(-1);
    gTruthArr.emplace_back(1);
    gTruthArr.emplace_back(0);

    int cnt = -1;
    int numErr = 0;
    for (std::string str : lhsStrArr) {
        ++cnt; 
        // Initialization
        Matrix<R> lhs;
        Matrix<R>::parseToMatrix(lhsStrArr[cnt], lhs);
        try {
            // Function call
            R det = lhs.getDeterminant();
            // Condition
            numErr += util::assertEq(det, gTruthArr[cnt], lhs, "matrix determinant");
        } catch (std::invalid_argument ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            ++numErr;
        }
    }
    if (cnt == -1) std::cout << "Warning: Did not test.\n";
    return numErr;
}

template<typename R>
int testNullspace() {
    std::vector<std::string> lhsStrArr;
    std::vector<int> lhsDimArr;
    std::vector<int> rankArr;

    lhsStrArr.emplace_back("1 4 0 2 3 6 0 6 7");
    lhsStrArr.emplace_back("0 0 1 0 1 0 1 0 0");
    lhsStrArr.emplace_back("1 0 2 0 1 0 2 0 4");
    lhsStrArr.emplace_back("0");
    lhsStrArr.emplace_back("1 1 0 4 3 0 2 3 2 5 5 4 0 4 6");
    lhsStrArr.emplace_back("3 3 0 0 2 3 6 2 -4 0 -12 0 0 3 0 2");
    lhsStrArr.emplace_back("3 1 -1 3 0 1 1 0");
    
    lhsDimArr.emplace_back(3);
    lhsDimArr.emplace_back(3);
    lhsDimArr.emplace_back(3);
    lhsDimArr.emplace_back(1);
    lhsDimArr.emplace_back(3);
    lhsDimArr.emplace_back(4);
    lhsDimArr.emplace_back(2);

    rankArr.emplace_back(0);
    rankArr.emplace_back(0);
    rankArr.emplace_back(1);
    rankArr.emplace_back(1);
    rankArr.emplace_back(2);
    rankArr.emplace_back(1);
    rankArr.emplace_back(2);
    
    int cnt = -1;
    int numErr = 0;
    for (auto x : lhsStrArr) {
        ++cnt;
        // Initialization
        Matrix<R> lhs;
        Matrix<R>::parseToMatrix(lhsStrArr[cnt], lhs, lhsDimArr[cnt]);

        try {
            // Function call
            Matrix<R> result = lhs.nullspace();

            Matrix<R> zero = lhs * result;
            int rank = result.getRank();

            // Assertion
            numErr += util::assertEq(zero, 0, lhs, "nullspace (validity))");
            numErr += util::assertEq(rank, rankArr[cnt], lhs, "nullspace (linear independence)");
        } catch (std::invalid_argument const& ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            ++numErr;
        }
    }
    if (cnt == -1) std::cout << "Warning: Did not test.\n";
    return numErr;
}


template<typename R>
int testCharPoly() {
    std::vector<std::string> lhsStrArr;
    std::vector<std::string> gTruthStrArr;

    lhsStrArr.emplace_back("-4 -4 8 -4 2 4 8 4 -4");
    lhsStrArr.emplace_back("1 -4 5 -4 1 4 8 4 -3");

    gTruthStrArr.emplace_back("224 -96 6 1");
    gTruthStrArr.emplace_back("219 -77 1 1");

    int cnt = -1;
    int numErr = 0;
    for (auto x : lhsStrArr) {
        ++cnt;
        // Initialization
        Matrix<R> lhs;
        Polynomial<R> gTruth;
        Matrix<R>::parseToMatrix(lhsStrArr[cnt], lhs);
        Polynomial<R>::parseToPolynomial(gTruthStrArr[cnt], gTruth);

        try {
            // Function call
            Polynomial<R> result = lhs.charPoly();
            // Assertion
            numErr += util::assertEq(result, gTruth, lhs, "characteristic polyomial");
        } catch (std::invalid_argument const& ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            ++numErr;
        }
    }
    if (cnt == -1) std::cout << "Warning: Did not test.\n";
    return numErr;
}

template<typename R>
int testComplete() {
    std::vector<std::string> lhsStrArr, rhsStrArr;
    std::vector<int> lhsDimArr;

    lhsStrArr.emplace_back("1 1");
    lhsStrArr.emplace_back("1 1 1 1 1 1 1 1 1");

    rhsStrArr.emplace_back("1 -4 5 -4 1 4");
    rhsStrArr.emplace_back("1 0 0 1 1 1 0 1 0");

    lhsDimArr.emplace_back(2);
    lhsDimArr.emplace_back(3);

    int cnt = -1;
    int numErr = 0;
    for (auto x : lhsStrArr) {
        ++cnt;
        // Initialization
        Matrix<R> lhs, rhs;
        Matrix<R>::parseToMatrix(lhsStrArr[cnt], lhs, lhsDimArr[cnt]);
        Matrix<R>::parseToMatrix(rhsStrArr[cnt], rhs, lhsDimArr[cnt]);

        try {
            int fullRank = rhs.getRank();
            // Function call
            Matrix<R> res = lhs.complete(rhs);
            int rankSum = lhs.getRank() + res.getM();
            Matrix<R> span = lhs.mergeHorizontal(res);
            int spanRank = span.getRank();

            // Assertion
            numErr += util::assertEq(rankSum, fullRank, lhs, "matrix completion (rank sum)");
            numErr += util::assertEq(spanRank, fullRank, lhs, "matrix completion (full rank)");
        } catch (std::invalid_argument const& ex) {
            std::cerr << ex.what() << " exception using " << typeid(R).name() << " coefficients.\n";
            ++numErr;
        }
    }
    if (cnt == -1) std::cout << "Warning: Did not test.\n";
    return numErr;
}
