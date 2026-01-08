#pragma once
#include "datastructure/matrix.hpp"
#include "datastructure/polynomial.hpp"

#include "util/helper.hpp"

#include <cassert>

namespace matOps {
	template <typename R>
	const auto rref(const Matrix<R>&, bool);

	template <typename R> const Matrix<R> inverse(const Matrix<R>&);
	template <typename R> const Polynomial<R> charPoly(const Matrix<R>&);
	template <typename R> const int rank(const Matrix<R>&);
	template <typename R> const R determinant(const Matrix<R>&);
	template <typename R> const Matrix<R> kernelBasis(const Matrix<R>&);
	template <typename R> const Matrix<R> completeBasis(const Matrix<R>&, const Matrix<R>&);
	template <typename R> const Matrix<R> minkSum(const Matrix<R>&, const Matrix<R>&);
	template <typename R> const Matrix<R> getVec(const Matrix<R>&, int);
}

/*
* @brief An implementation of the gaussian elimination algorithm, not requiring multiplicative inverses
* @param mat the Matrix to reduce
* @param reduced Returns echelon form if false, and reduced row echelon form if true
* 
* @return the echelon form, the rank and the determinant
*/
template <typename R>
const auto matOps::rref(const Matrix<R>& mat, bool reduced) {
	int n = mat.getN();
	int m = mat.getM();
	R* arr = util::allocate<R>(n * m);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			arr[m * i + j] = mat(i, j);
		}
	}
	R detNumer = 1;
	R detDenom = 1;
	int rank = 0;
	// The current number of reduced columns
	int numReduced = 0;
	// The row index of the current non-zero element
	int nzRowInd = -1;

	for (int j = 0; j < m; ++j) {
		nzRowInd = -1;
		// Find the first non-zero entry and its row-index s in column j
		for (int i = numReduced; i < n; ++i) {
			if (arr[m * i + j] != 0) {
				nzRowInd = i;
				break;
			}
		}
		// If there is no non zero entry, move to the next column
		if (nzRowInd == -1) {
			detDenom = 0;
			continue;
		}
		++rank;
		// Eradicate all other non-zero entries in column j
		for (int i = (reduced ? 0 : numReduced); i < n; ++i) {
			if (i != nzRowInd && arr[m * i + j] != 0) {
				R gcd = util::gcd(arr[m * nzRowInd + j], arr[m * i + j]);
				R nzRowFact = arr[m * i + j] / gcd;
				R elRowFact = arr[m * nzRowInd + j] / gcd;
				// Reduce row i and keep it as simple as possible
				for (int k = 0; k < m; ++k) {
					arr[m * i + k] *= elRowFact;
					arr[m * i + k] -= arr[m * nzRowInd + k] * nzRowFact;
				}
				detDenom *= elRowFact;
			}
			// Calculate the gcd of row i and factor it out
			R rowGcd = 0;
			for (int k = 0; k < m; ++k) {
				if (arr[m * i + k] != 0) {
					if (rowGcd == 0) {
						rowGcd = arr[m * i + k];
					}
					else {
						rowGcd = util::gcd(arr[m * i + k], rowGcd);
					}
				}
			}
			if (rowGcd != 0 && rowGcd != 1) {
				detNumer *= rowGcd;
				for (int k = 0; k < m; ++k) {
					arr[m * i + k] /= rowGcd;
				}
			}
		}
		// Switch rows nzRowInd and num
		if (nzRowInd != numReduced) {
			for (int k = 0; k < m; ++k) {
				std::swap(arr[m * nzRowInd + k], arr[m * numReduced + k]);
			}
			detNumer *= (-1);
		}
		++numReduced;
	}
	R det = 0;
	if (n == m && detDenom != 0) {
		for (int i = 0; i < n; ++i) {
			detNumer *= arr[m * i + i];
		}
		det = detNumer / detDenom;
	}
	struct result { Matrix<R> ref; int rank; R det; };
	result res{ Matrix<R>(arr, n, m), rank, det };
	util::deallocate(arr);
	return res;
}

/*
* @brief Calculate the inverse matrix or a scaled version of it in case of a ring
* Calculation is done by row reducing the matrix (mat, I) to (diag, mat^(-1)) where
* diag is a diagonal matrix
*/ 
template <typename R>
const Matrix<R> matOps::inverse(const Matrix<R>& mat) {
	int n = mat.getN();
	int m = mat.getM();
	if (n != m) {
		throw std::invalid_argument("Cannot find inverse of a non-square matrix.");
	}
	Matrix<R> I(1, n);
	// Get (mat, I)
	Matrix<R> help = minkSum(mat, I);
	// Reduce to (diag, mat^(-1))
	auto rrefRet = rref(help, true);
	Matrix<R> result = rrefRet.ref;
	// Find least common multiple of elements in the diagonal matrix diag
	R diagLcm = 1;
	for (int i = 0; i < n; ++i) {
		diagLcm *= result(i, i) / util::gcd(diagLcm, result(i, i));
	}
	// If diagLcm is zero, the rank of the matrix does not have full rank
	if (diagLcm == 0) {
		throw std::invalid_argument("Matrix is not invertible.");
	}
	// Scale mat^(-1) according to the lcm and the entry in the diagonal matrix in the corresponding row
	R* arr = util::allocate<R>(n * m);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			arr[i * m + j] = result(i, m + j) * diagLcm / result(i, i);
		}
	}
	result = Matrix<R>(arr, n, m);
	util::deallocate<R>(arr);
	//assert(mat * result == diagLcm);
	return result;
}

// Calculate the characteristic polynomial of a given matrix
template <typename R>
const Polynomial<R> matOps::charPoly(const Matrix<R>& mat) {
	int n = mat.getN();
	int m = mat.getM();
	if (n != m) {
		throw std::invalid_argument("Cannot calculate characteristic polynomial of non-quadratic matrix.");
	}
	Matrix<R> C = mat;
	R* arr = util::allocate<R>(n + 1);
	arr[n] = 1;
	for (int k = 1; k <= n; ++k) {
		if (k > 1) {
			C += Matrix<R>(arr[n - k + 1], n);
			C = mat * C;
		}
		arr[n - k] = 0;
		for (int j = 0; j < n; ++j) {
			arr[n - k] += C(j, j);
		}
		arr[n - k] /= (-k);
	}
	Polynomial<R> result(arr, n);
	util::deallocate(arr);
	return result;
}

// Return the rank
template <typename R>
const int matOps::rank(const Matrix<R>& mat) {
	return rref(mat, false).rank;
}

// Return the determinant
template <typename R>
const R matOps::determinant(const Matrix<R>& mat) {
	if (mat.getN() != mat.getM()) {
		throw std::invalid_argument("Cannot calculate determinant of non-square matrix.");
	}
	return rref(mat, false).det;
}

/*
* @brief Calculate a basis for the kernel of a given matrix
*/
template<typename R>
const Matrix<R> matOps::kernelBasis(const Matrix<R>& phi) {
	int n = phi.getN();
	int m = phi.getM();
	auto help = rref(phi, true);
	Matrix<R> phiRref = help.ref;
	if (help.rank == m) {
		throw std::invalid_argument("Cannot calculate basis, phi is injective.");
	}
	// Row indicates the current row in the rref, column counts the current number of vectors in the incomplete basis
	int row = 0;
	int column = 0;
	int dim = m - help.rank;
	// The indices of the first non-zero element in the rref
	int* indices = util::allocate<int>(help.rank);
	R* arr = util::allocate<R>(m * dim);

	for (int k = 0; k < m * dim; ++k) {
		arr[k] = 0;
	}

	for (int j = 0; j < m; ++j) {
		// if A[row][j] == 0, solve for variable x_j
		if (row >= n || phiRref(row, j) == 0) {
			// Calculate LCM of non-zero elements in column to avoid multiplicative inverses and unnecessarily big entries
			R colLcm(1);
			for (int k = 0; k < row; ++k) {
				if (phiRref(k, j) != 0) {
					colLcm *= phiRref(k, indices[k]) / util::gcd(colLcm, phiRref(k, indices[k]));
				}
			}
			// Add the vector corresponding to x_j to the result 
			arr[dim * j + column] = (colLcm == 0 ? 1 : colLcm);
			for (int k = 0; k < std::min(row + 1, n); ++k) {
				if (phiRref(k, j) != 0) {
					arr[dim * indices[k] + column] = phiRref(k, j) * colLcm / phiRref(k, indices[k]) * (-1);
				}
			}
			++column;
		} else {
			indices[row++] = j;
		}
	}
	Matrix<R> result(arr, m, dim);
	util::deallocate(arr);
	util::deallocate(indices);
	assert(phi * result == 0);
	return result;
}

/*
* Complete the column vectors of first to a basis of second. This is done by checking for each vector of second, if it is
* linearly independent to all vectors of this. If a linearly independent vector is found, add it to the basis and stop
* if the ranks of the two bases are equal
*/
template <typename R>
const Matrix<R> matOps::completeBasis(const Matrix<R>& first, const Matrix<R>& second) {
	int n = first.getN();
	int m = first.getM();
	if (first == 0) return second;
	// TODO: Check if span(this) is a subset of span(other)
	if (matOps::rank(first) == matOps::rank(second) || m == second.getM()) {
		return Matrix<R>(nullptr, n, 0);
	}
	// Check for every column vector in other, if it is linearly independent to the vectors in this basis.
	// If so, add it to the basis and continue until the ranks of the bases are equal.
	int newM = second.getM() - m;
	if (newM < 0) throw std::invalid_argument("Dimension of second larger than dimension of first.");
	R* arr = util::allocate<R>(n * second.getM());
	// arr = A^T, dimension m * n, dim - m row-vectors remain
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			arr[n * j + i] = first(i, j);
		}
	}
	int s = 0;
	int j = 0;
	do {
		// Add a vector of cspan as a row to the matrix
		int cRow = m + s + 1;
		for (int i = 0; i < n; ++i) {
			arr[n * (cRow - 1) + i] = second(i, j);
		}
		// Check if the rank has increased, else replace the vector by the next
		Matrix<R> help(arr, cRow, n);
		if (matOps::rank(help) > m + s) {
			++s;
		}
		++j;
	} while (s < newM);

	R* resArr = util::allocate<R>(n * newM);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < newM; ++j) {
			resArr[newM * i + j] = arr[n * (m + j) + i];
		}
	}
	util::deallocate(arr);
	// Reduce the columns
	for (int j = 0; j < newM; ++j) {
		R colGcd(0);
		for (int i = 0; i < n; ++i) {
			colGcd = (colGcd == 0 ? resArr[newM * i + j] : util::gcd(resArr[newM * i + j], colGcd));
		}
		if (colGcd != 0 && colGcd != 1) {
			for (int i = 0; i < n; ++i) {
				resArr[newM * i + j] /= colGcd;
			}
		}
	}
	Matrix<R> result(resArr, n, newM);
	util::deallocate(resArr);
	return result;
}

/*
* @brief Return the minkowski sum of the two matrices
* 
* A fancy way of returning (A, B) 
*/
template <typename R>
const Matrix<R> matOps::minkSum(const Matrix<R>& spanA, const Matrix<R>& spanB) {
	int n = spanA.getN();
	if (n != spanB.getN()) {
		throw std::invalid_argument("Cannot add vectors of different dimensions.");
	}
	int mA = spanA.getM();
	int mB = spanB.getM();
	int newM = mA + mB;
	R* arr = util::allocate<R>(n * newM);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < mA; ++j) {
			arr[newM * i + j] = spanA(i, j);
		}
		for (int j = 0; j < mB; ++j) {
			arr[newM * i + mA + j] = spanB(i, j);
		}
	}
	Matrix<R> result(arr, n, newM);
	util::deallocate(arr);
	return result;
}

// @brief return the j-th vector of mat
template<typename R>
const Matrix<R> matOps::getVec(const Matrix<R>& mat, const int j) {
	int n = mat.getN();
	R* arr = util::allocate<R>(n);
	for (int i = 0; i < n; ++i) {
		arr[i] = mat(i, j);
	}
	Matrix<R> result(arr, n, 1);
	util::deallocate<R>(arr);
	return result;
}