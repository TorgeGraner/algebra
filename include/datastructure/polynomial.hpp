#pragma once
#include <iostream>
#include <cassert>
#include "util/helper.hpp"

template <typename R> class Polynomial;

template <typename R> std::ostream& operator<<(std::ostream&, const Polynomial<R>&);
template <typename R> void swap(Polynomial<R>&, Polynomial<R>&);

template <typename R> 
class Polynomial {
private:
	R* coeffs = nullptr;
	int degree = -1;

	bool polyDiv(const Polynomial&, const Polynomial&, Polynomial&) const;	// Recursive polynomial division
	void normalize();								// Divide by gcd of coefficients

public:
	//-------------------------------------------------------------------------------------------------------------
	// Constructors
	//------------------------------------------//-----------------------------------------------------------------
	Polynomial() = default;						// Default constructor
	Polynomial(R*, const int);					// Standard constructor
	Polynomial(R, const int = 0);				// Monomial constructor
	Polynomial(int);							// Constant integer polynomial constructor
	Polynomial(const Polynomial&);				// Copy constructor
	Polynomial(Polynomial&&) noexcept;			// Move constructor

	~Polynomial() { util::deallocate(coeffs); }	// Destructor
	//------------------------------------------//-----------------------------------------------------------------
	// Operators
	//-------------------------------------------------------------------------------------------------------------
	Polynomial operator+=(const Polynomial&);
	Polynomial operator-=(const Polynomial&);
	Polynomial operator*=(const Polynomial&);
	Polynomial operator/=(const Polynomial&);
	Polynomial operator%=(const Polynomial&);

	friend Polynomial operator+(Polynomial lhs, const Polynomial& rhs) { return lhs += rhs; }
	friend Polynomial operator-(Polynomial lhs, const Polynomial& rhs) { return lhs -= rhs; }
	friend Polynomial operator*(Polynomial lhs, const Polynomial& rhs) { return lhs *= rhs; }
	friend Polynomial operator/(Polynomial lhs, const Polynomial& rhs) { return lhs /= rhs; }
	friend Polynomial operator%(Polynomial lhs, const Polynomial& rhs) { return lhs %= rhs; }

	Polynomial& operator=(const Polynomial&);		// Assignment operator
	Polynomial& operator=(Polynomial&&) noexcept;	// Move operator

	friend bool operator==(const Polynomial& lhs, const Polynomial& rhs) {
		if (lhs.degree != rhs.degree) return false;
		for (int k = 0; k <= lhs.degree; ++k)
			if (lhs[k] != rhs[k]) return false;
		return true;
	}
	const R& operator[](const int) const;			// Access operator

	friend bool operator!=(const Polynomial& lhs, const Polynomial& rhs) { return !(lhs == rhs); }
	friend std::ostream& operator<< <>(std::ostream&, const Polynomial&);
	//----------------------------------------------//-------------------------------------------------------------
	// Others
	//------------------------------------------//-----------------------------------------------------------------
	int getDegree() const { return degree; };	// Get the degree
	R map(const R&) const;						// Return value of this at specific value
	R coeffGcd() const;							// Return gcd of all coefficients
	Polynomial derivative() const;				// Return the algebraic derivative

	friend void swap<>(Polynomial&, Polynomial&);
};

//-------------------------------------------------------------------------------------------------------------
// Private
//-------------------------------------------------------------------------------------------------------------
/*
* Input: Polynomial q
* Output: Polynomial p such that *this = p * q + r with deg(r) < deg(*this)
* 
* By basic polynomial division, calculate a polynomial result, such that *this = divisor * result + rest.
* The rest in this will possess a degree less than result.
* 
* WARNING: Can, even in low dimensions, lead to an exponentially growing number of digits, overflowing standard int
*/
template<typename R>
bool Polynomial<R>::polyDiv(const Polynomial<R>& dividend, const Polynomial& divisor, Polynomial& res) const {
	int degF = dividend.getDegree();
	int degS = divisor.getDegree();
	int newDeg = degF - degS;

	if (degS < 0 || newDeg < 0) {
		return false;
	}
	if (degF < 0) {
		res = 0;
		return true;
	}

	R* resArr = util::allocate<R>(newDeg + 1);
	for (int i = 0; i < newDeg + 1; ++i) resArr[i] = 0;
	Polynomial rem = dividend;
	R remLead;
	R divisorLead = divisor[degS];

	while (newDeg >= 0) {
		// Eliminate the leading coefficient of the remaining polynomial 
		remLead = rem[degF];
		// If the lead of divisor does not divide the lead, the divisor polynomial does not divide this polynomial in this ring
		if (remLead % divisorLead != 0) {
			util::deallocate(resArr);
			return false;
		}
		resArr[newDeg] = remLead / divisorLead;
		// Enforce 0 = lead(remainder) * X^deg(remainder) - (lead(divisor) * X^deg(divisor)) * (lead(remainder) / lead(divisor)) * X^(deg(remainder) - deg(divisor))
		rem -= Polynomial(resArr[newDeg], newDeg) * divisor;
		// Should not occur and indicates faultily implemented ring class or integer overflow. Would result in a stack overflow
		assert(rem.getDegree() < degF);
		degF = rem.getDegree();
		newDeg = degF - degS;
	}
	res = Polynomial(resArr, degree - degS);
	util::deallocate(resArr);
	return true;
}

// Divide all coefficients by coefficient gcd
template<typename R>
void Polynomial<R>::normalize() {
	if (degree != -1) {
		R normalizer = coeffGcd();
		for (int k = 0; k <= degree; ++k) coeffs[k] /= normalizer;
	}
}

//-------------------------------------------------------------------------------------------------------------
// Constructors
//-------------------------------------------------------------------------------------------------------------
// Standard constructor
template <typename R>
Polynomial<R>::Polynomial(R* _coeffs, const int deg) : degree(deg) {
	if (degree != -1) {
		coeffs = util::allocate<R>(degree + 1);
		std::copy(_coeffs, _coeffs + degree + 1, coeffs);
	}
}

// Monomial constructor
template <typename R>
Polynomial<R>::Polynomial(R coeff, const int deg) : degree(deg) {
	if (coeff == 0) {
		degree = -1;
	} else {
		coeffs = util::allocate<R>(degree + 1);
		for (int k = 0; k < degree; ++k) {
			coeffs[k] = 0;
		}
		coeffs[degree] = coeff;
	}
}

// Constant integer polynomial constructor
template <typename R>
Polynomial<R>::Polynomial(int coeff) : degree(0) {
	if (coeff == 0) {
		degree = -1;
	} else {
		coeffs = util::allocate<R>(1);
		coeffs[0] = coeff;
	}
}

// Copy constructor
template <typename R>
Polynomial<R>::Polynomial(const Polynomial& orig) : degree(orig.degree) {
	if (degree != -1) {
		coeffs = util::allocate<R>(degree + 1);
		std::copy(orig.coeffs, orig.coeffs + degree + 1, coeffs);
	}
}

// Move constructor
template <typename R>
Polynomial<R>::Polynomial(Polynomial&& src) noexcept : Polynomial{} {
	swap(*this, src);
}

//-------------------------------------------------------------------------------------------------------------
// Operators
//-------------------------------------------------------------------------------------------------------------
template <typename R>
Polynomial<R> Polynomial<R>::operator+=(const Polynomial& addend) {
	int oDeg = addend.getDegree();
	int newDeg = std::max(degree, oDeg);
	// Calculate new degree
	if (degree == oDeg) {
		for (int k = degree; k >= 0; --k) {
			if ((*this)[k] + addend[k] != 0) break;
			--newDeg;
		}
	}
	R* arr = nullptr;
	if (newDeg != -1) {
		arr = util::allocate<R>(newDeg + 1);
		if (newDeg <= degree) {
			std::copy(coeffs, coeffs + newDeg + 1, arr);
			for (int k = 0; k <= std::min(oDeg, newDeg); ++k) {
				arr[k] += addend[k];
			}
		}
		else {
			// oDeg >= newDeg
			std::copy(addend.coeffs, addend.coeffs + newDeg + 1, arr);
			for (int k = 0; k < std::min(degree, newDeg) + 1; ++k) {
				arr[k] += (*this)[k];
			}
		}
	}
	util::deallocate(coeffs);
	degree = newDeg;
	coeffs = arr;
	return *this;
}

template <typename R>
Polynomial<R> Polynomial<R>::operator-=(const Polynomial& subtrahend) {
	return *this += subtrahend * (-1);
}

template <typename R>
Polynomial<R> Polynomial<R>::operator*=(const Polynomial& multiplicand) {
	int oDeg = multiplicand.getDegree();
	int newDeg = 0;

	if (degree == -1 || oDeg == -1) {
		newDeg = -1;
	} else {
		newDeg = degree + oDeg;
	}
	R* arr = nullptr;
	if (newDeg >= 0) {
		arr = util::allocate<R>(newDeg + 1);
		for (int i = 0; i < newDeg + 1; ++i) {
			arr[i] = 0;
			for (int j = 0; j <= i; ++j) {
				if (j < degree + 1 && i - j <= oDeg) {
					arr[i] += (*this)[j] * multiplicand[i - j];
				}
			}
		}
	}
	util::deallocate(coeffs);
	degree = newDeg;
	coeffs = arr;
	return *this;
}

template <typename R>
Polynomial<R> Polynomial<R>::operator/=(const Polynomial& divisor) {
	Polynomial q;
	if (*this == 0) return 0;
	if (!polyDiv(*this, divisor, q)) {
		std::cout << *this << " " << divisor << "\n";
		throw std::invalid_argument("Error in polynomial division.");
	}
	*this = q;
	return *this;
}

template <typename R>
Polynomial<R> Polynomial<R>::operator%=(const Polynomial& modulus) {
	if (modulus.getDegree() == -1) {
		throw std::invalid_argument("Cannot modulate by zero.");
	}
	if (*this == modulus) {
		*this = 0;
		return *this;
	}

	Polynomial<R> res;
	if (modulus.getDegree() == 0 || !polyDiv(*this, modulus, res)) {
		 // Polydiv cannot (or should not) be performed, enforce gcd(*this, modulus) == 1
		*this = (modulus == 1 ? 0 : 1);
	} else {
		// Correct but probably slows down the program considerably
		*this -= modulus * res;
	}
	return *this;
}

// Access operator for coefficients
template <typename R>
const R& Polynomial<R>::operator[](const int idx) const {
	if (idx < 0 || idx > degree) {
		throw std::invalid_argument("Index out of bounds for polynomial.");
	}
	return coeffs[idx];
}

// Assignment operator
template <typename R>
Polynomial<R>& Polynomial<R>::operator=(const Polynomial& rhs) {
	Polynomial tmp(rhs);
	swap(*this, tmp);
	return *this;
}

// Move operator
template <typename R>
Polynomial<R>& Polynomial<R>::operator=(Polynomial&& src) noexcept {
	Polynomial tmp(src);
	swap(*this, src);
	return *this;
}

//-------------------------------------------------------------------------------------------------------------
// Others
//-------------------------------------------------------------------------------------------------------------
// Return value at given value
template <typename R>
R Polynomial<R>::map(const R& x) const {
	R ret = 0;
	for (int i = 0; i <= degree; ++i) {
		ret *= x;
		ret += (*this)[degree - i];
	}
	return ret;
}

// Return algebraic derivative
template <typename R>
Polynomial<R> Polynomial<R>::derivative() const {
	if (degree < 1) return 0;
	R* arr = util::allocate<R>(degree);
	for (int k = 1; k <= degree; ++k) {
		arr[k - 1] = k * (*this)[k];
	}
	Polynomial<R> result(arr, degree - 1);
	util::deallocate(arr);
	return result;
}

// Return gcd of all coefficients
template <typename R>
R Polynomial<R>::coeffGcd() const {
	if (degree == -1) return 0;
	R res = coeffs[degree];
	for (int k = degree - 1; k >= 0; --k) {
		if (coeffs[k] != 0) {
			res = util::gcd(res, coeffs[k]);
		}
	}
	return res;
}

// Print operator
template <typename R>
std::ostream& operator<< <>(std::ostream& os, const Polynomial<R>& obj) {
	const int deg = obj.getDegree();
	// Zero polynomial
	if (deg == -1) {
		return os << 0;
	}
	if (deg == 0) {
		return os << obj[0];
	}
	bool first = true;
	for (int i = deg; i >= 0; --i) {
		if (obj[i] != 0) {
			if (!first) {
				os << "+";
			}
			if (i == 0 || obj[i] != 1) {
				os << obj[i];
			}
			if (i >= 1) {
				os << "X";
				if (i > 1) {
					os << "^" << i;
				}
			}
			first = false;
		}
	}
	return os;
}

template <typename R>
void swap(Polynomial<R>& lhs, Polynomial<R>& rhs) {
	std::swap(lhs.degree, rhs.degree);
	std::swap(lhs.coeffs, rhs.coeffs);
}

