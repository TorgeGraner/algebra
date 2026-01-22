#include <iostream>
#include "util/helper.hpp"
template <typename R> class Fractionalize;
template <typename R> std::ostream& operator<<(std::ostream&, const Fractionalize<R>&);
template <typename R>
/*
* @brief A datastructure implementing the field of fractions with coefficients in R
*/
class Fractionalize {
	private:
		R numerator;
		R denominator;

		void reduce();

	public:
		//-------------------------------------------------------------------------------------------------------------|
		// Constructors
		//-------------------------------------------------------------------------------------------------------------|
		Fractionalize() = default;
		Fractionalize(R);
		Fractionalize(int);
		Fractionalize(R, R);

		//-------------------------------------------------------------------------------------------------------------|
		// Operators
		//-------------------------------------------------------------------------------------------------------------|
		Fractionalize operator+=(const Fractionalize&);
		Fractionalize operator-=(const Fractionalize&);
		Fractionalize operator*=(const Fractionalize&);
		Fractionalize operator/=(const Fractionalize&);
		Fractionalize operator%=(const Fractionalize&);

		friend Fractionalize operator+(Fractionalize lhs, const Fractionalize& rhs) { return lhs += rhs; }
		friend Fractionalize operator-(Fractionalize lhs, const Fractionalize& rhs) { return lhs -= rhs; }
		friend Fractionalize operator*(Fractionalize lhs, const Fractionalize& rhs) { return lhs *= rhs; }
		friend Fractionalize operator/(Fractionalize lhs, const Fractionalize& rhs) { return lhs /= rhs; }
		friend Fractionalize operator%(Fractionalize lhs, const Fractionalize& rhs) { return lhs %= rhs; }

		bool operator==(const Fractionalize<R>&) const;
		bool operator!=(const Fractionalize<R>& rhs) const { return !(*this == rhs); }

		//-------------------------------------------------------------------------------------------------------------|
		// Getters
		//-------------------------------------------------------------------------------------------------------------|
		const R getNumerator() const { return numerator; }
		const R getDenominator() const { return denominator; }
};

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
	numerator *= rhs.denominator;
	denominator *= rhs.numerator;
	reduce(); 
	return *this;
}

template<typename R>
Fractionalize<R> Fractionalize<R>::operator%=(const Fractionalize& rhs) {
	numerator = 0;
	return *this;
}

template<typename R>
bool Fractionalize<R>::operator==(const Fractionalize& rhs) const {
	if (numerator == 0) return rhs.numerator == 0;
	return numerator == rhs.numerator && denominator == rhs.denominator;
}

template <typename R>
std::ostream& operator<< <>(std::ostream& os, const Fractionalize<R>& obj) {
	if (obj.getNumerator() == 0) return os << 0;
	if (obj.getDenominator() == 1) return os << obj.getNumerator();
	return os << obj.getNumerator() << "/" << obj.getDenominator();
}