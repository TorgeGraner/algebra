#pragma once

//-------------------------------------------------------------------------------------------------------------|
// Private
//---------------------------------------------------------------------|---------------------------------------|
// Polynomial division
//---------------------------------------------------------------------|

/*
* Standard polynomial division (iteratively, in place)
* Input: dividend, divisor, quotient and remainder such that
* divisor = dividend * quotient + remainder with deg(remainder) < deg(quotient)
* Output: Indicates if division was successful
*
* WARNING: Can, even in low dimensions, lead to integer overflow
*/
template<typename R>
bool Polynomial<R>::polyDiv(const Polynomial& divisor, Polynomial& quotientOut, Polynomial& remainderOut) const {
	int degF = degree;
	int degS = divisor.getDegree();

	int remDeg = degF;
	int resDeg = degF - degS;

	if (degS == -1 || resDeg < 0) {
		return false;
	}
	if (degF == -1) {
		quotientOut = 0;
		remainderOut = 0;
		return true;
	}

	// Coefficients of the result
	R* remArr = util::copy<R, R>(coeffs, degF + 1);
	R* resArr = util::zeroes<R>(resDeg + 1);

	R divLead = divisor(degS);

	for (int k = degF; k >= 0; --k) {
		R lead = remArr[k];
		if (lead != 0 && k < degS) break;
		--remDeg;
		if (lead == 0) continue;

		if (lead % divLead != 0) {
			util::deallocate(remArr);
			util::deallocate(resArr);
			return false;
		}
		// Calculate rem = rem - (rem_k / lead(divisor)) * X^(k - deg(divisor)) * divisor
		resArr[k - degS] = remArr[k] / divLead;
		for (int i = 0; i <= degS; ++i) {
			remArr[i + k - degS] -= resArr[k - degS] * divisor(i);
		}
	}
	quotientOut = Polynomial(resArr, resDeg);
	remainderOut = Polynomial(remArr, remDeg);
	return true;
}

//-------------------------------------------------------------------------------------------------------------|
// Public
//---------------------------------------------------------------------|---------------------------------------|
// Constructors
//---------------------------------------------------------------------|

// Monomial constructor
template <typename R>
Polynomial<R>::Polynomial(const R& coeff, int deg) : degree(deg) {
	if (coeff == 0) {
		degree = -1;
	} else {
		coeffs = util::zeroes<R>(degree + 1);
		coeffs[degree] = coeff;
	}
}

// Constant integer polynomial constructor
template <typename R>
Polynomial<R>::Polynomial(int coeff) : Polynomial<R>(R(coeff), 0) { }

// Copy constructor
template <typename R>
Polynomial<R>::Polynomial(const Polynomial& orig) : degree(orig.degree) {
	if (degree != -1) {
		coeffs = util::copy<R, R>(orig.coeffs, degree + 1);
	}
}

// Move constructor
template <typename R>
Polynomial<R>::Polynomial(Polynomial&& src) noexcept : Polynomial{} {
	swap(*this, src);
}

//---------------------------------------------------------------------|
// Operators (in place)
//---------------------------------------------------------------------|

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

//---------------------------------------------------------------------|
// Operators (out of place)
//---------------------------------------------------------------------|

template <typename R>
Polynomial<R> Polynomial<R>::operator+(const Polynomial& addend) const {
	int newDeg = std::max(degree, addend.getDegree());
	// Calculate new degree
	if (degree == addend.getDegree()) {
		while (newDeg >= 0 && (*this)(newDeg) + addend(newDeg) == 0) {
			--newDeg;
		}
	}
	if (newDeg == -1) return 0;

	R* arr = util::zeroes<R>(newDeg + 1);
	for (int i = 0; i <= std::min(newDeg, degree); ++i) arr[i] += (*this)(i);
	for (int i = 0; i <= std::min(newDeg, addend.getDegree()); ++i) arr[i] += addend(i);
	
	return Polynomial<R>(arr, newDeg);
}

template <typename R>
Polynomial<R> Polynomial<R>::operator-(const Polynomial& subtrahend) const {
	int newDeg = std::max(degree, subtrahend.getDegree());
	// Calculate new degree
	if (degree == subtrahend.getDegree()) {
		while (newDeg >= 0 && (*this)(newDeg) - subtrahend(newDeg) == 0) {
			--newDeg;
		}
	}
	if (newDeg == -1) return 0;

	R* arr = util::zeroes<R>(newDeg + 1);
	for (int i = 0; i <= std::min(newDeg, degree); ++i) arr[i] += (*this)(i);
	for (int i = 0; i <= std::min(newDeg, subtrahend.getDegree()); ++i) arr[i] -= subtrahend(i);
	
	return Polynomial<R>(arr, newDeg);
}

template <typename R>
Polynomial<R> Polynomial<R>::operator*(const Polynomial& multiplicand) const {
	if (degree == -1 || multiplicand.getDegree() == -1) return 0;

	int newDeg = degree + multiplicand.getDegree();
	R* arr = util::zeroes<R>(newDeg + 1);
	for (int i = 0; i <= degree; ++i) {
		for (int j = 0; j <= multiplicand.getDegree(); ++j) {
			arr[i + j] += (*this)(i) * multiplicand(j);
		}
	}
	return Polynomial<R>(arr, newDeg);
}

template <typename R>
Polynomial<R> Polynomial<R>::operator/(const Polynomial& divisor) const {
	if (divisor == 0) {
		throw std::invalid_argument("Division by zero.");
	}
	if (*this == 0) return 0;
	if (*this == divisor) return 1;

	// TODO: catch div not being a divisor
	Polynomial quotient, rem;
	if (!polyDiv(divisor, quotient, rem) || rem.getDegree() != -1) {
		throw std::invalid_argument("Cannot divide in this ring.");
	}
	return quotient;
}

template <typename R>
Polynomial<R> Polynomial<R>::operator%(const Polynomial& modulus) const {
	if (modulus == 0) {
		throw std::invalid_argument("Polynomial modulation by zero.");
	}
	if (*this == modulus) return 0;

	Polynomial<R> quotient, rem;
	if (modulus.getDegree() == 0 || !polyDiv(modulus, quotient, rem)) {
		// Polydiv cannot (or should not) be performed, assume gcd(*this, modulus) == 1
		return (modulus == 1 ? 0 : 1);
	} else {
		// Correct but probably slows down the program considerably
		return rem;
	}
}

// Polynomial comparison
template<typename R>
bool Polynomial<R>::operator==(const Polynomial& rhs) const {
	if (degree != rhs.degree) return false;
	for (int k = 0; k <= degree; ++k)
		if ((*this)(k) != rhs(k)) return false;
	return true;
}

// Access operator for coefficients
template <typename R>
const R& Polynomial<R>::operator()(int idx) const {
	if (idx < 0 || idx > degree) {
		throw std::invalid_argument("Index out of bounds.");
	}
	return coeffs[idx];
}

// Stream operator
template <typename R>
std::ostream& operator<< <>(std::ostream& os, const Polynomial<R>& obj) {
	const int deg = obj.getDegree();
	// Zero polynomial
	if (deg == -1) {
		return os << 0;
	}
	if (deg == 0) {
		return os << obj(0);
	}
	bool first = true;
	for (int i = deg; i >= 0; --i) {
		if (obj(i) != 0) {
			if (!first) {
				os << "+";
			}
			if (i == 0 || obj(i) != 1) {
				os << obj(i);
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

//---------------------------------------------------------------------|
// Operations
//---------------------------------------------------------------------|


// Return value at given value using horners method
template <typename R>
R Polynomial<R>::map(const R& x) const {
	R ret = 0;
	for (int i = 0; i <= degree; ++i) {
		ret *= x;
		ret += (*this)(degree - i);
	}
	return ret;
}

// Return algebraic derivative
template <typename R>
Polynomial<R> Polynomial<R>::derivative() const {
	if (degree < 1) return 0;
	R* arr = util::allocate<R>(degree);
	for (int k = 1; k <= degree; ++k) {
		arr[k - 1] = k * (*this)(k);
	}
	return Polynomial<R>(arr, degree - 1);
}

//---------------------------------------------------------------------|
// Others
//---------------------------------------------------------------------|

template <typename R>
void swap(Polynomial<R>& lhs, Polynomial<R>& rhs) {
	std::swap(lhs.degree, rhs.degree);
	std::swap(lhs.coeffs, rhs.coeffs);
}

//---------------------------------------------------------------------|
// Static functions
//---------------------------------------------------------------------|

// Parse string of coefficients into polynomial
template <typename R> 
bool Polynomial<R>::parseToPolynomial(const std::string& str, Polynomial<R>& out) {
	int num = 0;
	std::istringstream is(str);
    std::vector<int> data;
    
	while (is >> num) data.emplace_back(num);

	int numTotal = int(size(data));
    
    R* values = util::allocate<R>(numTotal);
	int cnt = 0;
    for (int value : data) values[cnt++] = value;

    out = Polynomial<R>(values, numTotal - 1);
	return true;
}
