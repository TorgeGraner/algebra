#pragma once
#include "datastructure/polynomial.hpp"
#include "datastructure/matrix.hpp"

#include "algorithm/matrix_operations.hpp"

#include "util/helper.hpp"

/*
* The algebraic behaviour of a root of a polynomial is completely determined by the polynomial. Even though only
* the polynomial might be known, the minimal polynomial of the sum and product of such roots can be calculated, yielding
* a field. 
* WARNING: Contains a conceptual flaw, which needs to get fixed and which leads to x - x not necessarily being zero.
*/
template <typename R>
class Root {
private:
	Polynomial<R> minPoly;
	int id = 0;

public:
	//-------------------------------------------------------------------------------------------------------------|
	// Constructors
	//-------------------------------------------------------------------------------------------------------------|
	Root() = default;
	Root(const Polynomial<R>&, const int);
	Root(const R&);
	Root(const int);
	//-------------------------------------------------------------------------------------------------------------|
	// Operators
	//-------------------------------------------------------------------------------------------------------------|
	Root operator+=(const Root&);
	Root operator-=(const Root&);
	Root operator*=(const Root&);
	Root operator/=(const Root&);

	bool operator==(const Root&) const;
	bool operator!=(const Root&) const;
	//-------------------------------------------------------------------------------------------------------------|
	// Others
	//-------------------------------------------------------------------------------------------------------------|
	const Polynomial<R> getMinPoly() const { return minPoly; }
	const int getId() const { return id; }

	friend Root operator+(Root lhs, const Root& rhs) { return lhs += rhs; }
	friend Root operator-(Root lhs, const Root& rhs) { return lhs -= rhs; }
	friend Root operator*(Root lhs, const Root& rhs) { return lhs *= rhs; }
	friend Root operator/(Root lhs, const Root& rhs) { return lhs /= rhs; }

	friend std::ostream& operator<<(std::ostream& str, const Root& obj) {
		return str << "{" << obj.minPoly << " = 0, x_" << obj.id << "}";
	}
};

template <typename R>
Root<R>::Root(const Polynomial<R>& minPoly, const int id) : minPoly(minPoly), id(id) {
	if (minPoly.getDegree() < 1) {
		throw std::invalid_argument("Constant Polynomial does not possess any zeros.");
	}
	if (id < 0 || id > minPoly.getDegree()) {
		throw std::invalid_argument("ID invalid for Root of Polynomial.");
	}
}

template <typename R>
Root<R>::Root(const R& val) {
	minPoly = Polynomial<R>(1, 1) - val;
	id = 0;
}

template <typename R>
Root<R>::Root(const int val) {
	minPoly = Polynomial<R>(1, 1) - R(val);
	id = 0;
}

/*
* The shared minimal polynomial of ALL possible sums of roots of this and addend
*/
template <typename R>
Root<R> Root<R>::operator+=(const Root& addend) {
	Polynomial<R> oPoly = addend.getMinPoly();

	int degA = minPoly.getDegree();
	int degB = oPoly.getDegree();
	int nDeg = degA * degB;

	// The polynomials that replace the leading exponents. Note X^degA = mPolA and Y^degB = mPolB
	Polynomial<R> mPolA(minPoly[degA], degA);
	Polynomial<R> mPolB(oPoly[degB], degB);
	mPolA -= minPoly;
	mPolB -= oPoly;

	Matrix<R> span(1, nDeg, 1);
	// The coefficient of X^j * Y^k is stored in arr[j * degB + k]
	R* arr = util::allocate<R>(nDeg);
	for (int i = 1; i < nDeg + 1; ++i) {
		for (int k = 0; k < nDeg; ++k) arr[k] = 0;
		// Get next row by multiplying Elements of last row by (x + y) 
		for (int j = 0; j < degA; ++j) {
			for (int k = 0; k < degB; ++k) {
				// Get X^j * Y^k
				R entry = span(j * degB + k, i - 1); 
				if (entry == 0) continue;
				// Generate the powers X^{j+1} * Y^k and X^j * Y^{k+1}
				// Generate X^{j+1} * Y^k 
				if (j + 1 == degA) { 
					// Case X^{j+1} = mPolA, so mPolA * Y^k with powers X^l * Y^k
					for (int l = 0; l <= mPolA.getDegree(); ++l) {
						arr[l * degB + k] += entry * mPolA[l];
					}
				} else { 
					// Case X^{j+1} * Y^k
					arr[(j + 1) * degB + k] += entry;
				}
				// Generate X^j * Y^{k+1}
				if (k + 1 == degB) {
					// Case Y^{k+1} = mPolB, so X^j * mPolB with powers X^j * Y^l
					for (int l = 0; l <= mPolB.getDegree(); ++l) {
						arr[j * degB + l] += entry * mPolB[l];
					}
				} else { 
					// Case X^j * Y^{k+1}
					arr[j * degB + (k + 1)] += entry;
				}
			}
		}
		span = matOps::minkSum(span, Matrix<R>(arr, nDeg, 1));

		auto help = matOps::rref(span, false);
		// Check if the generated powers are already linearly dependent
		if (help.rank < i) {
			span = help.ref;
			break;
		}
	}
	util::deallocate(arr);
    // A linear dependence was found, the minimal polynomial calculates from the kernel of the span of the powers
	span = matOps::kernelBasis(span);
	arr = util::allocate<R>(span.getN());
	for (int i = 0; i < span.getN(); ++i) {
		arr[i] = span(i, 0);
	}
	minPoly = Polynomial<R>(arr, span.getN() - 1);
	id = 0;																									// TODO: Fix id
	util::deallocate(arr);
	return *this;
}

template <typename R>
Root<R> Root<R>::operator-=(const Root& subtrahend) {
	return *this += subtrahend * (-1);
}

template <typename R>
Root<R> Root<R>::operator*=(const Root& multiplicand) {
	Polynomial<R> oPoly = multiplicand.getMinPoly();

	int degA = minPoly.getDegree();
	int degB = oPoly.getDegree();
	int nDeg = degA * degB;

	// The polynomials that replace the leading exponents. Note X^degA = mPolA and Y^degB = mPolB
	Polynomial<R> mPolA(minPoly[degA], degA);
	Polynomial<R> mPolB(oPoly[degB], degB);
	mPolA -= minPoly;
	mPolB -= oPoly;

	Matrix<R> span(1, nDeg, 1);
	// The coefficient of X^j * Y^k is stored in arr[j * degB + k]
	R* arr = util::allocate<R>(nDeg);
	for (int i = 1; i < nDeg + 1; ++i) {
		for (int k = 0; k < nDeg; ++k) arr[k] = 0;
		// Get next row by multiplying Elements of last row by XY 
		for (int j = 0; j < degA; ++j) {
			for (int k = 0; k < degB; ++k) {
				// Get X^j * Y^k
				R entry = span(j * degB + k, i - 1); 
				if (entry == 0) continue;
				// Generate the power X^{j+1} * Y^{k+1}
				if (j + 1 == degA) { 
					if (k + 1 == degB) {
						// Case mPolA * mPolB with powers X^a * Y^b
						for (int a = 0; a <= mPolA.getDegree(); ++a)
							for (int b = 0; b <= mPolB.getDegree(); ++b)
								arr[a * degB + b] += entry * mPolA[a] * mPolB[b];
					} else {
						// Case mPolA * Y^{k+1} with powers X^a * Y^{k+1}
						for (int a = 0; a <= mPolA.getDegree(); ++a)
							arr[a * degB + (k + 1)] += entry * mPolA[a];
					}
				} else {
					if (k + 1 == degB) {
						// Case X^{j+1} * mPolB
						for (int b = 0; b <= mPolB.getDegree(); ++b)
							arr[(j + 1) * degB + b] += entry * mPolB[b];
					} else {
						// Case X^{j+1} * Y^{k+1}
						arr[(j + 1) * degB + (k + 1)] += entry;
					}
				}
			}
		}
		span = matOps::minkSum(span, Matrix<R>(arr, nDeg, 1));
		auto help = matOps::rref(span, false);
		// Check if the generated powers are already linearly dependent
		if (help.rank < i) {
			span = help.ref;
			break;
		}
	}
	util::deallocate(arr);
    // A linear dependence was found, the minimal polynomial calculates from the kernel of the span of the powers
	span = matOps::kernelBasis(span);
	arr = util::allocate<R>(span.getN());
	for (int i = 0; i < span.getN(); ++i) {
		arr[i] = span(i, 0);
	}
	minPoly = Polynomial<R>(arr, span.getN() - 1);
	id = 0;																									// TODO: Fix id
	util::deallocate(arr);
	return *this;
}

/*
* Do the same as with multiplication
*/
template <typename R>
Root<R> Root<R>::operator/=(const Root& dividend) {
	Polynomial<R> oPoly = dividend.getMinPoly();

	int degA = minPoly.getDegree();
	int degB = oPoly.getDegree();
	int nDeg = degA * degB;
	
	Polynomial<R> x(1, 1);

	// Find the inverse of Y
	// Divide by biggest possible power of Y
	int l = 0;
	while (oPoly[l] == 0) ++l;
	R lCoeff = oPoly[l];
	Polynomial<R> yInv = minPoly / Polynomial<R>(-1, l);
	// Write a_l/Y in terms of powers of Y, i.e. a_l/Y = -sum_{k=l+1}^na_kY^{k-l-1}
	yInv -= yInv[0];
	yInv /= x;

	// The polynomials that replace the leading exponents. Note X^degA = mPolA
	Polynomial<R> mPolA(minPoly[degA], degA);
	mPolA -= minPoly;

	Matrix<R> span(1, nDeg, 1);
	R* arr = util::allocate<R>(nDeg);
	for (int i = 1; i < nDeg + 1; ++i) {
		// Allocate array of zeros. The coefficient of X^j * Y^k is stored in arr[j * degB + k]
		for (int k = 0; k < nDeg; ++k) arr[k] = 0;
		// Get next row by multiplying Elements of last row by X * Y^{-1}
		for (int j = 0; j < degA; ++j) {
			for (int k = 0; k < degB; ++k) {
				// Get X^j * Y^k
				R entry = span(j * degB + k, i - 1); 
				if (entry == 0) continue;
				// Generate the power X^{j+1} * Y^{k-1}
				if (j + 1 == degA) { 
					if (k == 0) {
						// Case mPolA * yInv with powers X^a * Y^b
						for (int a = 0; a <= mPolA.getDegree(); ++a)
							for (int b = 0; b <= yInv.getDegree(); ++b)
								arr[a * degB + b] += entry * mPolA[a] * yInv[b];
					} else {
						// Case mPolA * Y^{k-1} with powers X^a * Y^{k-1}
						for (int a = 0; a <= mPolA.getDegree(); ++a)
							arr[a * degB + (k - 1)] += lCoeff * entry * mPolA[a];
					}
				} else {
					if (k == 0) {
						// Case X^{j+1} * yInv
						for (int b = 0; b <= yInv.getDegree(); ++b)
							arr[(j + 1) * degB + b] += entry * yInv[b];
					} else {
						// Case X^{j+1} * Y^{k-1}
						arr[(j + 1) * degB + (k - 1)] += lCoeff * entry;
					}
				}
			}
		}
		// If lCoeff != 1, the preceding rows need to be scaled accordingly, such that the linear dependence gives the correct polynomial
		span *= lCoeff;
		span = matOps::minkSum(span, Matrix<R>(arr, nDeg, 1));

		auto help = matOps::rref(span, false);
		// Check if the generated powers are already linearly dependent
		if (help.rank < i) {
			span = help.ref;
			break;
		}
	}
	util::deallocate(arr);
    // A linear dependence was found, the minimal polynomial calculates from the kernel of the span of the powers
	span = matOps::kernelBasis(span);
	arr = util::allocate<R>(span.getN());
	for (int i = 0; i < span.getN(); ++i) {
		arr[i] = span(i, 0);
	}
	minPoly = Polynomial<R>(arr, span.getN() - 1);
	id = 0;																									// TODO: Fix id
	util::deallocate(arr);
	return *this;
}

// TODO: Test for correctness
template <typename R>
bool Root<R>::operator==(const Root& other) const {
	return util::gcd(minPoly, other.minPoly) != 1;
}

template <typename R>
bool Root<R>::operator!=(const Root& other) const {
	return !(*this == other);
}