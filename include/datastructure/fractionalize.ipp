
/**
 * @brief Reduces the numerator and denominator of the fraction by their greatest common divisor.
 * 
 * @tparam R The type of underlying ring
 */
template<typename R>
void Fractionalize<R>::reduce() {
	R g = util::gcd(numerator, denominator);
	if (g != 0) {
		numerator /= g;
		denominator /= g;
	}
}

template<typename R>
Fractionalize<R>::Fractionalize(int num) : numerator(num), denominator(1) {};

template<typename R> 
Fractionalize<R>::Fractionalize(R num) : numerator(num), denominator(1) {};

template<typename R> 
Fractionalize<R>::Fractionalize(R num, R denom) : numerator(num), denominator(denom) { reduce(); };

template<typename R>
Fractionalize<R> Fractionalize<R>::operator+=(const Fractionalize& rhs) {
	numerator *= rhs.denominator;
	numerator += denominator * rhs.numerator;
	denominator *= rhs.denominator;
	reduce();
	return *this;
}

template<typename R>
Fractionalize<R> Fractionalize<R>::operator-=(const Fractionalize& rhs) {
	numerator *= rhs.denominator;
	numerator -= denominator * rhs.numerator;
	denominator *= rhs.denominator;
	reduce();
	return *this;
}

template<typename R>
Fractionalize<R> Fractionalize<R>::operator*=(const Fractionalize& rhs) {
	numerator *= rhs.numerator;
	denominator *= rhs.denominator;
	reduce();
	return *this;
}

template<typename R>
Fractionalize<R> Fractionalize<R>::operator/=(const Fractionalize& rhs) {
	if (rhs == 0) throw std::invalid_argument("Division by zero.");
	numerator *= rhs.denominator;
	denominator *= rhs.numerator;
	reduce(); 
	return *this;
}

template<typename R>
Fractionalize<R> Fractionalize<R>::operator%=(const Fractionalize& rhs) {
	if (rhs == 0) throw std::invalid_argument("Modulation by zero.");
	numerator = 0;
	return *this;
}

template<typename R>
bool Fractionalize<R>::operator==(const Fractionalize& rhs) const {
	if (numerator == 0) return rhs.numerator == 0;
	return numerator == rhs.numerator && denominator == rhs.denominator;
}

/**
 * @brief Output operator for fractions
 * 
 * @tparam R The type of underlying ring
 * @param os the output stream
 * @param obj the fraction object to be printed
 * 
 * @return the output stream
 */
template <typename R>
std::ostream& operator<< <>(std::ostream& os, const Fractionalize<R>& obj) {
	if (obj.getNumerator() == 0) return os << 0;
	if (obj.getDenominator() == 1) return os << obj.getNumerator();
	return os << obj.getNumerator() << "/" << obj.getDenominator();
}