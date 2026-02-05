#pragma once

//-------------------------------------------------------------------------------------------------------------|
// Private
//---------------------------------------------------------------------|---------------------------------------|
// Polynomial division
//---------------------------------------------------------------------|

/**
* @brief Standard polynomial division (iteratively)

* @tparam R The type of the underlying ring
* @param[in] rhs: The polynomial to divide by
* @param[out] quotientOut: The resulting quotient
* @param[out] remainderOut: The resulting remainder

* @warning Can, even in low dimensions, lead to integer overflow

* @return True if the division was successful, false if not 
*/
template<typename R>
bool Polynomial<R>::polyDiv(const Polynomial& rhs, Polynomial& quotientOut, Polynomial& remainderOut) const {
	int degF = degree;
	int degS = rhs.getDegree();

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

	R divLead = rhs(degS);

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
		// Calculate rem = rem - (rem_k / lead(rhs)) * X^(k - deg(rhs)) * rhs
		resArr[k - degS] = remArr[k] / divLead;
		for (int i = 0; i <= degS; ++i) {
			remArr[i + k - degS] -= resArr[k - degS] * rhs(i);
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

template <typename R>
Polynomial<R>::Polynomial(const R& coeff, int deg) : degree(deg) {
	if (coeff == 0) {
		degree = -1;
	} else {
		coeffs = util::zeroes<R>(degree + 1);
		coeffs[degree] = coeff;
	}
}

template <typename R>
Polynomial<R>::Polynomial(int coeff) : Polynomial<R>(R(coeff), 0) {}

template <typename R>
Polynomial<R>::Polynomial(const Polynomial& orig) : degree(orig.degree) {
	if (degree != -1) {
		coeffs = util::copy<R, R>(orig.coeffs, degree + 1);
	}
}

template <typename R>
Polynomial<R>::Polynomial(Polynomial&& src) noexcept : Polynomial{} {
	swap(*this, src);
}

//---------------------------------------------------------------------|
// Operators (in place)
//---------------------------------------------------------------------|

template <typename R>
Polynomial<R>& Polynomial<R>::operator=(const Polynomial& rhs) {
	Polynomial tmp(rhs);
	swap(*this, tmp);
	return *this;
}

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
Polynomial<R> Polynomial<R>::operator+(const Polynomial& rhs) const {
	int newDeg = std::max(degree, rhs.getDegree());
	// Calculate new degree
	if (degree == rhs.getDegree()) {
		while (newDeg >= 0 && (*this)(newDeg) + rhs(newDeg) == 0) {
			--newDeg;
		}
	}
	if (newDeg == -1) return 0;

	R* arr = util::zeroes<R>(newDeg + 1);
	for (int i = 0; i <= std::min(newDeg, degree); ++i) arr[i] += (*this)(i);
	for (int i = 0; i <= std::min(newDeg, rhs.getDegree()); ++i) arr[i] += rhs(i);
	
	return Polynomial<R>(arr, newDeg);
}

template <typename R>
Polynomial<R> Polynomial<R>::operator-(const Polynomial& rhs) const {
	int newDeg = std::max(degree, rhs.getDegree());
	// Calculate new degree
	if (degree == rhs.getDegree()) {
		while (newDeg >= 0 && (*this)(newDeg) - rhs(newDeg) == 0) {
			--newDeg;
		}
	}
	if (newDeg == -1) return 0;

	R* arr = util::zeroes<R>(newDeg + 1);
	for (int i = 0; i <= std::min(newDeg, degree); ++i) arr[i] += (*this)(i);
	for (int i = 0; i <= std::min(newDeg, rhs.getDegree()); ++i) arr[i] -= rhs(i);
	
	return Polynomial<R>(arr, newDeg);
}

template <typename R>
Polynomial<R> Polynomial<R>::operator*(const Polynomial& rhs) const {
	if (degree == -1 || rhs.getDegree() == -1) return 0;

	int newDeg = degree + rhs.getDegree();
	R* arr = util::zeroes<R>(newDeg + 1);
	for (int i = 0; i <= degree; ++i) {
		for (int j = 0; j <= rhs.getDegree(); ++j) {
			arr[i + j] += (*this)(i) * rhs(j);
		}
	}
	return Polynomial<R>(arr, newDeg);
}

/**
 * @brief Polynomial division
 * 
 * @tparam R The type of the underlying ring
 * @param rhs The polynomial to divide by
 * 
 * @throws std::invalid_argument if the division is not possible
 * 
 * @return The quotient of the two polynomials
 */
template <typename R>
Polynomial<R> Polynomial<R>::operator/(const Polynomial& rhs) const {
	if (rhs == 0) {
		throw std::invalid_argument("Division by zero.");
	}
	if (*this == 0) return 0;
	if (*this == rhs) return 1;

	Polynomial quotient, rem;
	if (!polyDiv(rhs, quotient, rem) || rem.getDegree() != -1) {
		throw std::invalid_argument("Cannot divide in this ring.");
	}
	return quotient;
}

/**
 * @brief Polynomial modulation
 * 
 * @tparam R The type of the underlying ring
 * @param rhs The polynomial to modulate by
 * 
 * @throws std::invalid_argument if the modulus is zero
 * 
 * @return The remainder of the division of this by rhs, or zero if polynomial division can (or should) not be performed 
 */
template <typename R>
Polynomial<R> Polynomial<R>::operator%(const Polynomial& rhs) const {
	if (rhs == 0) {
		throw std::invalid_argument("Polynomial modulation by zero.");
	}
	if (*this == rhs) return 0;

	Polynomial<R> quotient, rem;
	if (rhs.getDegree() == 0 || !polyDiv(rhs, quotient, rem)) {
		// Polydiv cannot (or should not) be performed, assume gcd(*this, rhs) == 1
		return (rhs == 1 ? 0 : 1);
	} else {
		return rem;
	}
}

template<typename R>
bool Polynomial<R>::operator==(const Polynomial& rhs) const {
	if (degree != rhs.degree) return false;
	for (int k = 0; k <= degree; ++k)
		if ((*this)(k) != rhs(k)) return false;
	return true;
}

/**
 * @brief Access operator for coefficients
 * 
 * @tparam R The type of the underlying ring
 * @param idx The index of the coefficient to access
 * 
 * @throws std::invalid_argument if the index is out of bounds
 * 
 * @return The coefficient at the given index
 */
template <typename R>
const R& Polynomial<R>::operator()(int idx) const {
	if (idx < 0 || idx > degree) {
		throw std::invalid_argument("Index out of bounds.");
	}
	return coeffs[idx];
}

/**
 * @brief Output operator for polynomials
 * 
 * @tparam R The type of the underlying ring
 * @param os The output stream
 * @param obj The polynomial object to be printed
 * 
 * @return The output stream
 */
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

/**
 * @brief Evaluates the polynomial at a given value using Horner's method
 * 
 * @tparam R The type of the underlying ring
 * @param x The value to evaluate the polynomial at
 * 
 * @return The value of the polynomial at x
 */
template <typename R>
R Polynomial<R>::map(const R& x) const {
	R ret = 0;
	for (int i = 0; i <= degree; ++i) {
		ret *= x;
		ret += (*this)(degree - i);
	}
	return ret;
}

/**
 * @brief Returns the algebraic derivative of the polynomial
 * 
 * @tparam R The type of the underlying ring
 * 
 * @return The derivative of the polynomial
 */
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

/**
 * @brief A parse function for polynomials, which reads coefficients from a string and constructs a polynomial from them. The coefficients are expected to be in the form of "a_n a_(n-1) ... a_0", where a_i is the coefficient of X^i.
 * 
 * @tparam R The type of the underlying ring
 * @param str The string to be parsed
 * @param out The polynomial to be constructed
 * 
 * @return true if the parsing was successful, false otherwise
 */
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
