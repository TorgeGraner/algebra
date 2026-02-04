#include <ostream>

#include "util/helper.hpp"

template <typename R> class Fractionalize;
template <typename R> std::ostream& operator<<(std::ostream&, const Fractionalize<R>&);

/*
* @brief A datastructure implementing the field of fractions with coefficients in R
*/
template <typename R>
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
		static int characteristic() { return R::characteristic(); }
};

#include "fractionalize.tpp"