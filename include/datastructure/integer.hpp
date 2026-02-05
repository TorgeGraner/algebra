#pragma once
#include <ostream>

class Integer;
std::ostream& operator<<(std::ostream&, const Integer&);

/**
* @brief A datastructure implementing the ring of integers
*/
class Integer {
private:
	int64_t value = 0;	///< The integer value

public:
	//-------------------------------------------------------------------------------------------------------------|
	// Constructors
	//-------------------------------------------------------------------------------------------------------------|

	Integer() = default;
	Integer(int64_t val) : value(val) {};	///< Constructor from integer

	//-------------------------------------------------------------------------------------------------------------|
	// Operators
	//-------------------------------------------------------------------------------------------------------------|

	Integer operator+=(const Integer& rhs) { value += rhs.value; return *this; }	///< In place addition
	Integer operator-=(const Integer& rhs) { value -= rhs.value; return *this; }	///< In place subtraction
	Integer operator*=(const Integer& rhs) { value *= rhs.value; return *this; }	///< In place multiplication
	Integer operator/=(const Integer& rhs) { 
		if (rhs == 0) throw std::invalid_argument("Cannot divide by zero.");
		value /= rhs.value; 
		return *this; 
	}														///< In place division
	Integer operator%=(const Integer& rhs) { 
		if (rhs == 0) throw std::invalid_argument("Cannot divide by zero.");
		value %= rhs.value; 
		return *this; 
	}														///< In place modulation

	friend Integer operator+(Integer lhs, const Integer& rhs) { return lhs += rhs; }	///< Out of place addition
	friend Integer operator-(Integer lhs, const Integer& rhs) { return lhs -= rhs; }	///< Out of place subtraction
	friend Integer operator*(Integer lhs, const Integer& rhs) { return lhs *= rhs; }	///< Out of place multiplication
	friend Integer operator/(Integer lhs, const Integer& rhs) { return lhs /= rhs; }	///< Out of place division
	friend Integer operator%(Integer lhs, const Integer& rhs) { return lhs %= rhs; }	///< Out of place modulation

	bool operator==(const Integer& rhs) const { return value == rhs.value; }	///< Equality operator
	bool operator!=(const Integer& rhs) const { return value != rhs.value; }	///< Inequality operator

	//-------------------------------------------------------------------------------------------------------------|
	// Getters
	//-------------------------------------------------------------------------------------------------------------|

	int getValue() const { return int(value); }	///< Get the integer value
	static int characteristic() { return 0; }	///< Get the characteristic
};

/**
 * @brief Stream operator
 * 
 * @param os The output stream
 * @param obj The object to stream
 * @return The output stream
 */
std::ostream& operator<<(std::ostream& os, const Integer& obj) { 
	return os << obj.getValue();
}