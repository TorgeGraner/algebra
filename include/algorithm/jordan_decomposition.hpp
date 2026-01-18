#pragma once
#include <vector>

#include "datastructure/matrix.hpp"
#include "datastructure/polynomial.hpp"

#include "algorithm/matrix_operations.hpp"

template <typename R> bool getFactor(const Polynomial<R>&, Polynomial<R>&);
template <typename R> Matrix<R> decompose(Matrix<R>&);


/*
* @brief Calculate a factor of a given polynomial
* @param toFactor the polynomial to factor
* @param factor A pointer to a polynomial, to which the result is written
* @return true if a factor was found, false otherwise
* 
* TODO: Improve everything
*/
template <typename R>
bool getFactor(const Polynomial<R>& toFactor, Polynomial<R>& factor) {
	if (toFactor.getDegree() == -1) return false;
	Polynomial<R> X(1, 1);
	// If the constant term is zero, X is a factor
	if (toFactor[0] == 0) {
		factor = X;
		return true;
	}
	// If toFactor is of the Form aX + b return X + b/a 
	if (toFactor.getDegree() == 1) {
		factor = X + Polynomial<R>(toFactor[0] / toFactor[1]);					// TODO: Possible error, rather return aX + b
		return true;
	}
	// If the degree if of toFactor is 2, return its roots by the abc-formula
	if (toFactor.getDegree() == 2) {
		const R& a = toFactor[2];
		const R& b = toFactor[1];
		const R& c = toFactor[0];

		if (b * b - a * c * 4 == 0) {
			factor = X + Polynomial<R>(b / (a * 2));							// TODO: Possible error
			return true;
		}
	}
	// Brute force check for p(n) = 0 for some integers n
	for (int tVal = -20; tVal <= 20; ++tVal) {
		R val(tVal);
		if (toFactor.map(val) == 0) {
			factor = X - Polynomial<R>(val);
			return true;
		}
	}
	return false;
}

/*
* Return the (possibly partial) Jordan normal form
* TODO: Return Frobenius normal form for factors of degree two
*/
template <typename R>
Matrix<R> decompose(Matrix<R>& A) {
	int n = A.getN();
	Polynomial<R> cPoly = matOps::charPoly(A);
	// Calculate factorization of the characteristic polynomial 
	Matrix<R> jordanBase = Matrix<R>(nullptr, n, 0);
	std::vector<R> eigenvalues;
	Polynomial<R> currFactor;

	while (getFactor(cPoly, currFactor)) {
		R eig = currFactor[0] * (-1);
		// Calculate the algebraic multiplicity of the current eigenvalue
		int alg_mult = 0;
		while (cPoly.map(eig) == 0) {
			cPoly /= currFactor;
			++alg_mult;
		}
		std::cout << "Found eigenvalue " << eig << " of algebraic multiplicity " << alg_mult << std::endl;
		std::vector<Matrix<R>> complements;
		// Set psi := A - lambda * I
		Matrix<R> psi = A - Matrix<R>(eig, n);
		Matrix<R> mult = Matrix<R>(1, n);
		Matrix<R> lastKernel = Matrix<R>(nullptr, n, 0);
		// Calculate the complements ker psi^(k+1) / ker psi^k, containing possible base vectors for the primary components
		for (int i = 0; i < alg_mult; ++i) {
			Matrix<R> help = matOps::kernelBasis(mult *= psi);
			complements.push_back(matOps::completeBasis(lastKernel, help));
			lastKernel = help;
		}
		for (int i = (int) complements.size() - 1; i >= 0; --i) {
			// Add all jordan chains of remaining vectors
			while (complements[i].getM() != 0) {
				Matrix<R> vect = matOps::getVec(complements[i], complements[i].getM() - 1);
				// Add the jordan chain corresponding to vect
				for (int k = 0; k < i + 1; ++k) {
					jordanBase = matOps::minkSum(vect, jordanBase);
					// Remove psi^k(v) from the corresponding complement, such that no vector in its span is chosen later
					complements[i - k] = matOps::completeBasis(vect, complements[i - k]);
					vect = psi * vect;
				}
			}
		}
	}
	return matOps::minkSum(jordanBase, matOps::completeBasis(jordanBase, Matrix<R>(1, n)));
}
