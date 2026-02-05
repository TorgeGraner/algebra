#pragma once
#include <ostream>
#include <stdexcept>

template <typename int p> class ModP;
template <typename int p> std::ostream& operator<<(std::ostream&, const ModP<p>&);

/**
 * @brief A class representing the integers modulo p
 * 
 * @note Will be a field iff p is prime
 * 
 * @tparam p The modulus
 */
template <typename int p>
class ModP {
private:
	int value = 0; 	///< The integer value 

	inline int mod(int x) const {
		while (x < 0) x += p;
		return x %= p;
	}				///< Modified modulation, ensuring non-negative results

public:
	//-------------------------------------------------------------------------------------------------------------|
	// Constructors
	//-------------------------------------------------------------------------------------------------------------|

	ModP() = default;
	ModP(int val) : value(mod(val)) {};		///< Constructor from integer

	//-------------------------------------------------------------------------------------------------------------|
	// Operators
	//-------------------------------------------------------------------------------------------------------------|

	ModP operator+=(const ModP& rhs) { value = mod(value + rhs.value); return *this; }	///< In place addition
	ModP operator-=(const ModP& rhs) { value = mod(value - rhs.value); return *this; }	///< In place subtraction
	ModP operator*=(const ModP& rhs) { value = mod(value * rhs.value); return *this; }	///< In place multiplication
	ModP operator%=(const ModP& rhs) { 
		if (rhs == 0) throw std::invalid_argument("Cannot modulate by zero.");
		value = 0; 					   
		return *this; 
	}	///< In place modulation
	ModP operator/=(const ModP& rhs) {
        if (rhs == 0) throw std::invalid_argument("Cannot divide by zero.");
        for (int i = 0; i < p; ++i) {
            if (i * rhs.value % p == value) {
                value = i;
                return *this;
            }
        }
		// Can only occur iff p is not prime and rhs a divisor of p
        throw std::invalid_argument("Element does not divide this.");
	}																			///< In place division

	friend ModP operator+(ModP lhs, const ModP& rhs) { return lhs += rhs; }		///< Out of place addition
	friend ModP operator-(ModP lhs, const ModP& rhs) { return lhs -= rhs; }		///< Out of place subtraction	
	friend ModP operator*(ModP lhs, const ModP& rhs) { return lhs *= rhs; }		///< Out of place multiplication
	friend ModP operator/(ModP lhs, const ModP& rhs) { return lhs /= rhs; }		///< Out of place division
	friend ModP operator%(ModP lhs, const ModP& rhs) { return lhs %= rhs; }		///< Out of place modulation

	bool operator==(const ModP& rhs) const { return value == rhs.value; }		///< Equality operator
	bool operator!=(const ModP& rhs) const { return !(*this == rhs); }			///< Inequality operator

	//-------------------------------------------------------------------------------------------------------------|
	// Getters
	//-------------------------------------------------------------------------------------------------------------|

	int getValue() const { return value; }		///< Get the integer value
	static int characteristic() { return p; }	///< Get the characteristic
};

/**
 * @brief Output stream operator
 * 
 * @tparam p The modulus
 * @param os The output stream
 * @param obj The object to stream
 * @return The output stream
 */
template<typename int p>
std::ostream& operator<<(std::ostream& os, const ModP<p>& obj) { 
	return os << "[" << obj.getValue() << "]";
}