#pragma once
#include <vector>

#include "datastructure/matrix.hpp"
#include "datastructure/polynomial.hpp"

template <typename R> bool getRoot(const Polynomial<R>&, Polynomial<R>&);
template <typename R> Matrix<R> jordanDecomposition(Matrix<R>&);

/**
 * @brief Get a root of a given polyomial.
 * 
 * @tparam R The type of the underlying ring
 * @param[in] toFactor The polynomial to get a root of
 * @param[out] root The root of the polynomial, if found
 * 
 * @return true if a root was found, false otherwise
 */
template <typename R>
bool getRoot(const Polynomial<R>& toFactor, R& root) {
	if (toFactor.getDegree() < 1) return false;
	// If the constant term is zero, X is a factor
	if (toFactor(0) == 0) {
		root = 0;
		return true;
	}
	for (int val = -20; val <= 20; ++val) {
		if (toFactor.map(val) == 0) {
			root = val;
			return true;
		}
	}
	// If toFactor is of the Form aX + b return X + b/a 
	if (toFactor.getDegree() == 1) {
		root = toFactor(0) / toFactor(1);					// TODO: Possible error, rather return aX + b
		return true;
	}
	// If the degree if of toFactor is 2, return its roots by the abc-formula
	if (toFactor.getDegree() == 2) {
		R a = toFactor(2);
		R b = toFactor(1);
		R c = toFactor(0);

		if (b * b - a * c * 4 == 0) {
			root = b / (a * 2);							// TODO: Possible error
			return true;
		}
	}
	return false;
}

/**
 * @brief Calculate a jordan basis of a given matrix
 * 
 * @tparam R The type of the underlying euclidean ring
 * @param A The matrix to calculate the jordan basis of
 * 
 * @return A matrix whose columns form a jordan basis of A
 */
template <typename R>
Matrix<R> jordanDecomposition(Matrix<R>& A) {
	int n = A.getN();
	Polynomial<R> X = Polynomial<R>(1, 1);
	Polynomial<R> cPoly = A.charPoly();

	Matrix<R> jordanBase = Matrix<R>(R(0), n, 0);
	R eig;
	while (getRoot(cPoly, eig)) {
		Polynomial<R> currFactor = X - eig;
		// Calculate the algebraic multiplicity of the current eigenvalue
		int algMult = 0;
		while (cPoly.map(eig) == 0) {
			cPoly /= currFactor;
			++algMult;
		}
		std::vector<Matrix<R>> kernels;
		// Set psi := A - lambda * I
		Matrix<R> psi = A - Matrix<R>(eig, n, n);
		// psi^k, initially psi^0 = I
		Matrix<R> psiPower = Matrix<R>(R(1), n, n);
		// ker(psi^k), initially ker(psi^0) = 0
		int lastKernelDimension = -1;
		// Calculate the complements ker psi^(k+1) / ker psi^k, containing possible base vectors for the primary components
		for (int i = 0; i <= algMult; ++i) {
			Matrix<R> currKernel = psiPower.nullspace();
			if (lastKernelDimension == currKernel.getM()) break;
			kernels.insert(kernels.begin(), currKernel);
			psiPower = psiPower * psi;
			lastKernelDimension = currKernel.getM();
		}
		int numBlocks = int(kernels.size()) - 1;
		// psi(ker psi^k) subseteq ker psi^{k - 1}
		Matrix<R> imageSum = Matrix<R>(R(0), n, 0);
		for (int k = 0; k < numBlocks; ++k) {
			kernels[k] = imageSum.mergeHorizontal(kernels[k + 1]).complete(kernels[k]);
			imageSum = psi * imageSum.mergeHorizontal(kernels[k]);
		}
		for(int k = 0; k < numBlocks; ++k) {
			// All remaining vectors generate unique jordan chains
			for (int j = 0; j < kernels[k].getM(); ++j) {
				Matrix<R> vect = kernels[k].getColumns(j, j);
				for (int p = 0; p < numBlocks - k; ++p) { 
					jordanBase = vect.mergeHorizontal(jordanBase);
					vect = psi * vect;
				}
			}
		}
	}
	if(jordanBase.getM() != n) {
		jordanBase = jordanBase.mergeHorizontal(jordanBase.complete(Matrix<R>(R(1), n, n)));
	}
	return jordanBase;
}
