#pragma once
#include <iostream>

class Integer;
std::ostream& operator<<(std::ostream&, const Integer&);
/*
* @brief A datastructure implementing the ring of integers, as a wrapper class of the native type int
* This also fixes the problem associated with int, which drops the fractional part in a division (i.e. int x = 5/2 
* where x now has the value 2, resulting in faulty logic)
*/
class Integer {
private:
	int value = 0;

public:
	//-------------------------------------------------------------------------------------------------------------|
	// Constructors
	//-------------------------------------------------------------------------------------------------------------|
	Integer() = default;
	Integer(int val) : value(val) {};

	//-------------------------------------------------------------------------------------------------------------|
	// Operators
	//-------------------------------------------------------------------------------------------------------------|
	Integer operator+=(const Integer& rhs) { value += rhs.value; return *this; }
	Integer operator-=(const Integer& rhs) { value -= rhs.value; return *this; }
	Integer operator*=(const Integer& rhs) { value *= rhs.value; return *this; }
	Integer operator/=(const Integer& rhs) { 
	if (value % rhs.value != 0) throw std::invalid_argument("Cannot divide by an integer that is not a divisor.");
		value /= rhs.value; 
		return *this; 
	}
	Integer operator%=(const Integer& rhs) { value %= rhs.value; return *this; }

	friend Integer operator+(Integer lhs, const Integer& rhs) { return lhs += rhs; }
	friend Integer operator-(Integer lhs, const Integer& rhs) { return lhs -= rhs; }
	friend Integer operator*(Integer lhs, const Integer& rhs) { return lhs *= rhs; }
	friend Integer operator/(Integer lhs, const Integer& rhs) { return lhs /= rhs; }
	friend Integer operator%(Integer lhs, const Integer& rhs) { return lhs %= rhs; }

	bool operator==(const Integer& rhs) const { return value == rhs.value; }
	bool operator!=(const Integer& rhs) const { return value != rhs.value; }

	//-------------------------------------------------------------------------------------------------------------|
	// Getters
	//-------------------------------------------------------------------------------------------------------------|
	int getValue() const { return value; }
};

std::ostream& operator<<(std::ostream& os, const Integer& obj) { 
	return os << obj.getValue();
}