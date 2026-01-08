#pragma once
#include <iostream>

class Integer;
std::ostream& operator<<(std::ostream&, const Integer&);

class Integer {
private:
	int value = 0;

public:
	//-------------------------------------------------------------------------------------------------------------|
	// Constructor
	//-------------------------------------------------------------------------------------------------------------|
	Integer() = default;	// Default constructor
	Integer(int val) : value(val) {};

	//-------------------------------------------------------------------------------------------------------------|
	// Operators
	//-------------------------------------------------------------------------------------------------------------|
	Integer operator+=(const Integer& rhs) { value += rhs.value; return *this; }
	Integer operator-=(const Integer& rhs) { value -= rhs.value; return *this; }
	Integer operator*=(const Integer& rhs) { value *= rhs.value; return *this; }
	Integer operator/=(const Integer& rhs) { 
		//if (value %= rhs.value != 0) throw std::invalid_argument("Cannot divide by an integer that is not a divisor.");
		value /= rhs.value; 
		return *this; 
	}
	Integer operator%=(const Integer& rhs) { value %= rhs.value; return *this; }

	bool operator==(const Integer& rhs) const { return value == rhs.value; }
	bool operator!=(const Integer& rhs) const { return value != rhs.value; }

	friend Integer operator+(Integer lhs, const Integer& rhs) { return lhs += rhs; }
	friend Integer operator-(Integer lhs, const Integer& rhs) { return lhs -= rhs; }
	friend Integer operator*(Integer lhs, const Integer& rhs) { return lhs *= rhs; }
	friend Integer operator/(Integer lhs, const Integer& rhs) { return lhs /= rhs; }
	friend Integer operator%(Integer lhs, const Integer& rhs) { return lhs %= rhs; }

	int getValue() const { return value; }
};

std::ostream& operator<<(std::ostream& os, const Integer& obj) { 
	return os << obj.getValue();
}