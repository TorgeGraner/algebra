#pragma once
#include <ostream>
#include <stdexcept>

template <typename int epsilon> class EpsFloat;
template <typename int epsilon> std::ostream& operator<<(std::ostream&, const EpsFloat<epsilon>&);

/**
 * @brief A class representing a floating point number, with equality defined up to a certain precision
 * 
 * @tparam epsilon The required precision, i.e. two EpsFloats a and b are considered equal iff |a - b| < 1/epsilon
 */
template <typename int epsilon>
class EpsFloat {
    private:
        double value = 0;	///< The double value of this EpsFloat

    public:

	//-------------------------------------------------------------------------------------------------------------|
	// Constructors
	//-------------------------------------------------------------------------------------------------------------|

	EpsFloat() = default;
	EpsFloat(double val) : value(val) {};	///< Constructor from double

	//-------------------------------------------------------------------------------------------------------------|
	// Operators
	//-------------------------------------------------------------------------------------------------------------|

	EpsFloat operator+=(const EpsFloat& rhs) { value += rhs.value; return *this; }	///< In place addition
	EpsFloat operator-=(const EpsFloat& rhs) { value -= rhs.value; return *this; }	///< In place subtraction
	EpsFloat operator*=(const EpsFloat& rhs) { value *= rhs.value; return *this; }	///< In place multiplication
	EpsFloat operator/=(const EpsFloat& rhs) { 
		if (rhs == 0) throw std::invalid_argument("Cannot divide by zero.");
		value /= rhs.value; 
		return *this; 
	}																				///< In place division
	EpsFloat operator%=(const EpsFloat& rhs) { 
		if (rhs == 0) throw std::invalid_argument("Cannot modulate by zero.");
		value = 0;          
		return *this; 
	}																				///< In place modulation

	friend EpsFloat operator+(EpsFloat lhs, const EpsFloat& rhs) { return lhs += rhs; }
	friend EpsFloat operator-(EpsFloat lhs, const EpsFloat& rhs) { return lhs -= rhs; }
	friend EpsFloat operator*(EpsFloat lhs, const EpsFloat& rhs) { return lhs *= rhs; }
	friend EpsFloat operator/(EpsFloat lhs, const EpsFloat& rhs) { return lhs /= rhs; }
	friend EpsFloat operator%(EpsFloat lhs, const EpsFloat& rhs) { return lhs %= rhs; }

	bool operator==(const EpsFloat& rhs) const { return std::abs(value - rhs.value) < (1 / double(epsilon)); }
	bool operator!=(const EpsFloat& rhs) const { return !(*this == rhs); }

	//-------------------------------------------------------------------------------------------------------------|
	// Getters
	//-------------------------------------------------------------------------------------------------------------|

	double getValue() const { return value; }
	static int characteristic() { return 0; }
};

template<typename int p>
std::ostream& operator<<(std::ostream& os, const EpsFloat<p>& obj) { 
    if (obj == 0) return os << 0;
    else return os << obj.getValue();
}