#pragma once
#include "datastructure/polynomial.hpp"
#include "datastructure/matrix.hpp"

#include "algorithm/matrix_operations.hpp"

#include "util/helper.hpp"

// The algebraic behaviour of a root of a polynomial is completely determined by the polynomial. Even though only
// the polynomial might be known, the minimal polynomial of the sum and product of such roots can be calculated, yielding
// a field. 
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

template <typename R>
Root<R> Root<R>::operator+=(const Root& addend) {
	Polynomial<R> oPoly = addend.getMinPoly();
	Polynomial<R> x(1, 1);
	Polynomial<R> xPol = 1, yPol = 1;

	int degA = minPoly.getDegree();
	int degB = oPoly.getDegree();
	int nDeg = degA * degB;

	// The polynomials that replace the leading exponents. Note X^degA = mPolA and Y^degB = mPolB
	Polynomial<R> mPolA(minPoly[degA], degA);
	Polynomial<R> mPolB(oPoly[degB], degB);
	mPolA -= minPoly;
	mPolB -= oPoly;

	Matrix<R> span(nullptr, nDeg, 0);
	R* arr = nullptr;
	for (int i = 0; i < nDeg + 1; ++i) {
		// Allocate vector of zeros
		arr = util::allocate<R>(nDeg);
		for (int k = 0; k < nDeg; ++k) arr[k] = 0;
        
        if (i == 0) { 
            arr[0] = 1;
        } else {    
            for (int j = 0; j < degA; ++j) {
                for (int k = 0; k < degB; ++k) {
                    R entry = span(j * degB + k, i - 1); // Get X^j * Y^k in span[j * degB + k] -> X^{j+1} * Y^k and X^j * Y^{k+1}
                    if (entry != 0) {
                        if (j + 1 == degA) { // -> X^{j+1} = mPolA, so mPolA * Y^k und X^l * Y^k
                            for (int l = 0; l <= mPolA.getDegree(); ++l) 
                                arr[l * degB + k] += entry * mPolA[l];
                        } else { // -> X^{j+1} * Y^k
                            arr[(j + 1) * degB + k] += entry;
                        }

                        if (k + 1 == degB) { // -> Y^{k+1} = mPolB, so X^j * mPolB und X^j * Y^l
                            for (int l = 0; l <= mPolB.getDegree(); ++l) 
                                arr[j * degB + l] += entry * mPolB[l];
                        } else { // -> X^j * Y^{k+1}
                            arr[j * degB + (k + 1)] += entry;
                        }
                    }
                }
            }
        }
        // Add a new column
		span = matOps::minkSum(span, Matrix<R>(arr, nDeg, 1));
		util::deallocate(arr);

		auto help = matOps::rref(span, false);
		// Check if the generated powers are already linearly dependent
		if (help.rank < i) {
			span = help.ref;
			break;
		}
	}
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
Root<R> Root<R>::operator*=(const Root& multiplicand) { // Inefficient
	Polynomial<R> oPoly = multiplicand.getMinPoly();
	Polynomial<R> x(1, 1);
	Polynomial<R> xPol = 1, yPol = 1;

	int degA = minPoly.getDegree();
	int degB = oPoly.getDegree();
	int nDeg = degA * degB;

	// Polynomials that replace the leading exponents
	Polynomial<R> mPolA(minPoly[degA], degA);
	Polynomial<R> mPolB(oPoly[degB], degB);
	mPolA -= minPoly;
	mPolB -= oPoly;

	Matrix<R> span(nullptr, nDeg, 0);
	R* arr = nullptr;
	for (int i = 0; i < nDeg + 1; ++i) {
		// Allocate vector of zeros
		arr = util::allocate<R>(nDeg);
		for (int k = 0; k < nDeg; ++k) arr[k] = 0;
		// Input the coefficient of the power x^j * y^k in row j * degB + k.
		for (int j = 0; j <= xPol.getDegree(); ++j) {
			for (int k = 0; k <= yPol.getDegree(); ++k) {
				arr[j * degB + k] = xPol[j] * yPol[k];
			}
		}
        // Add a new column
		span = matOps::minkSum(span, Matrix<R>(arr, nDeg, 1));
		util::deallocate(arr);
		auto help = matOps::rref(span, false);
		// Check if the generated powers are already linearly dependent
		if (help.rank == i) {
			span = help.ref;
			break;
		}
		// Generate the next power (xy)^{i + 1}
		xPol *= x;
		yPol *= x;
		R lead;
        // Replace highest powers corresponding to the minimal polynomials
		if (xPol.getDegree() == degA) {
			lead = xPol[degA];
			xPol -= Polynomial<R>(lead, degA);
			xPol += lead * mPolA;
		}
		if (yPol.getDegree() == degB) {
			lead = yPol[degB];
			yPol -= Polynomial<R>(lead, degB);
			yPol += lead * mPolB;
		}
	}
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
Root<R> Root<R>::operator/=(const Root& dividend) {
	Polynomial<R> oPoly = dividend.getMinPoly();
	Polynomial<R> x(1, 1);
	Polynomial<R> pol = 1;

	int degA = minPoly.getDegree();
	int degB = oPoly.getDegree();
	int nDeg = degA * degB;

	Polynomial<R> mPol(oPoly[degB], degB);
	mPol -= oPoly;

	Polynomial<R>* yPols = util::allocate<Polynomial<R>>(nDeg + 1); // to free
	R* arr = util::allocate<R>(nDeg * (nDeg + 1));
	for (int k = 0; k < nDeg * (nDeg + 1); ++k) arr[k] = 0;

	for (int i = 0; i < nDeg + 1; ++i) {
		yPols[i] = pol;
		pol *= x;
		R lead;
		if (pol.getDegree() == degB) {
			lead = pol[degB];
			pol -= Polynomial<R>(lead, degB);
			pol += lead * mPol;
		}
	}
	pol = 1;
	mPol = Polynomial<R>(minPoly[degA], degA);
	mPol -= minPoly;
	for (int i = 0; i < nDeg + 1; ++i) {
		Polynomial<R> yPol = yPols[nDeg - i];
		for (int j = 0; j <= pol.getDegree(); ++j) {
			for (int k = 0; k <= yPol.getDegree(); ++k) {
				arr[(j * degB + k) * (nDeg + 1) + i] = pol[j] * yPol[k];
			}
		}
		pol *= x;
		R lead;
		if (pol.getDegree() == degA) {
			lead = pol[degA];
			pol -= Polynomial<R>(lead, degA);
			pol += lead * mPol;
		}
	}

	Matrix<R> help(arr, nDeg, nDeg + 1);
	util::deallocate(arr);
	help = matOps::kernelBasis(help);
	arr = util::allocate<R>(nDeg + 1);
	for (int i = 0; i <= nDeg; ++i) {
		arr[i] = help(i, 0);
	}
	minPoly = Polynomial<R>(arr, nDeg);
	id = 0;																									// TODO: Fix id
	util::deallocate(yPols);
	util::deallocate(arr);
	return *this;
}

// TODO: Test for correctness
template <typename R>
bool Root<R>::operator==(const Root& other) const {
	return (minPoly == other.minPoly()) && (id == other.id());
}

template <typename R>
bool Root<R>::operator!=(const Root& other) const {
	return !(*this == other);
}