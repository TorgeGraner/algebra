#include <ostream>

#include "util/helper.hpp"

template <typename R> class Fractionalize;
template <typename R> std::ostream& operator<<(std::ostream&, const Fractionalize<R>&);

/**
* @brief A class representing a fraction of two elements 
* @tparam The underlying integral domain
*/
template <typename R>
class Fractionalize {
	private:
		R numerator;	///< The numerator of this fraction
		R denominator;	///< The denominator of this fraction

		void reduce();

	public:
		//-------------------------------------------------------------------------------------------------------------|
		// Constructors
		//-------------------------------------------------------------------------------------------------------------|
		
		Fractionalize() = default;
		Fractionalize(R);			///< Constructs a fraction with the given numerator
		Fractionalize(int);			///< Constructs a fraction with the given numerator
		Fractionalize(R, R);		///< Constructor from numerator and denominator

		//-------------------------------------------------------------------------------------------------------------|
		// Operators
		//-------------------------------------------------------------------------------------------------------------|

		Fractionalize operator+=(const Fractionalize&);	///< In place fraction addition
		Fractionalize operator-=(const Fractionalize&);	///< In place fraction subtraction
		Fractionalize operator*=(const Fractionalize&);	///< In place fraction multiplication
		Fractionalize operator/=(const Fractionalize&);	///< In place fraction division
		Fractionalize operator%=(const Fractionalize&);	///< In place fraction modulation

		friend Fractionalize operator+(Fractionalize lhs, const Fractionalize& rhs) { return lhs += rhs; }	// Out of place fraction addition
		friend Fractionalize operator-(Fractionalize lhs, const Fractionalize& rhs) { return lhs -= rhs; } 	// Out of place fraction subtraction
		friend Fractionalize operator*(Fractionalize lhs, const Fractionalize& rhs) { return lhs *= rhs; }	// Out of place fraction multiplication
		friend Fractionalize operator/(Fractionalize lhs, const Fractionalize& rhs) { return lhs /= rhs; }	// Out of place fraction division
		friend Fractionalize operator%(Fractionalize lhs, const Fractionalize& rhs) { return lhs %= rhs; }	// Out of place fraction modulation

		bool operator==(const Fractionalize<R>&) const;									///< Fraction equality
		bool operator!=(const Fractionalize<R>& rhs) const { return !(*this == rhs); }	///< Fraction inequality

		//-------------------------------------------------------------------------------------------------------------|
		// Getters
		//-------------------------------------------------------------------------------------------------------------|

		const R getNumerator() const { return numerator; }			///< Get the numerator
		const R getDenominator() const { return denominator; }		///< Get the denominator
		static int characteristic() { return R::characteristic(); }	///< Get the characteristic of the underlying integral domain
};

#include "fractionalize.ipp"