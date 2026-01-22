#pragma once
#include <vector>

#include "datastructure/matrix.hpp"
#include "datastructure/polynomial.hpp"

#include "algorithm/matrix_operations.hpp"

template <typename R> bool getRoot(const Polynomial<R>&, Polynomial<R>&);
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
bool getRoot(const Polynomial<R>& toFactor, R& root) {
	if (toFactor.getDegree() < 1) return false;
	Polynomial<R> X(1, 1);
	// If the constant term is zero, X is a factor
	if (toFactor[0] == 0) {
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
		root = toFactor[0] / toFactor[1];					// TODO: Possible error, rather return aX + b
		return true;
	}
	// If the degree if of toFactor is 2, return its roots by the abc-formula
	if (toFactor.getDegree() == 2) {
		const R& a = toFactor[2];
		const R& b = toFactor[1];
		const R& c = toFactor[0];

		if (b * b - a * c * 4 == 0) {
			root = b / (a * 2);							// TODO: Possible error
			return true;
		}
	}
	// Brute force check for p(n) = 0 for some integers n
	return false;
}

/*
* Return the (possibly partial) Jordan normal form
* TODO: Return Frobenius normal form for factors of degree two
*/
template <typename R>
Matrix<R> decompose(Matrix<R>& A) {
	int n = A.getN();

	Polynomial<R> X(1, 1);
	Polynomial<R> cPoly = matOps::charPolyNaive(A);
	std::vector<R> eigenvalues;

	Matrix<R> jordanBase = Matrix<R>(nullptr, n, 0);
	R eig;
	while (getRoot(cPoly, eig)) {
		Polynomial<R> currFactor = X - eig;
		// Calculate the algebraic multiplicity of the current eigenvalue
		int alg_mult = 0;
		while (cPoly.map(eig) == 0) {
			cPoly /= currFactor;
			++alg_mult;
		}
		//std::cout << "Found eigenvalue " << eig << " of algebraic multiplicity " << alg_mult << "\n";
		std::vector<Matrix<R>> kernels;
		// Set psi := A - lambda * I
		Matrix<R> psi = A - Matrix<R>(eig, n, n);
		// psi^k, initially psi^0 = I
		Matrix<R> psiPower = Matrix<R>(1, n, n);
		// ker(psi^k), initially ker(psi^0) = 0
		Matrix<R> lastKernel = Matrix<R>(nullptr, n, 0);
		// Calculate the complements ker psi^(k+1) / ker psi^k, containing possible base vectors for the primary components
		for (int i = 0; i < alg_mult; ++i) {
			psiPower *= psi;
			Matrix<R> currKernel = matOps::kernelBasis(psiPower);
			// TODO: Validate
			if (lastKernel.getM() == currKernel.getM()) break;
			kernels.insert(kernels.begin(), currKernel);
			lastKernel = currKernel;
		}
		int numBlocks = int(kernels.size());
		// psi(ker psi^k) subseteq ker psi^{k - 1}
		Matrix<R> imageSum(nullptr, n, 0);
		for (int k = 0; k < numBlocks - 1; ++k) {
			kernels[k] = matOps::completeBasis(matOps::minkSum(imageSum, kernels[k + 1]), kernels[k]);
			imageSum = psi * matOps::minkSum(imageSum, kernels[k]);
		}
		kernels[numBlocks - 1] = matOps::completeBasis(imageSum, kernels[numBlocks - 1]);

		for(int k = 0; k < numBlocks; ++k) {
			// All remaining vectors generate unique jordan chains
			for (int j = 0; j < kernels[k].getM(); ++j) {
				Matrix<R> vect = matOps::getVec(kernels[k], j);
				for (int p = 0; p < numBlocks - k; ++p) {
					jordanBase = matOps::minkSum(vect, jordanBase);
					vect = psi * vect;
				}
			}
		}
	}
	return matOps::minkSum(jordanBase, matOps::completeBasis(jordanBase, Matrix<R>(1, n, n)));
}
