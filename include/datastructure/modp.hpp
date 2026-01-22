#pragma once
#include <iostream>

template <typename int p> class ModP;
template <typename int p> std::ostream& operator<<(std::ostream&, const ModP<p>&);

/*
* @brief A datastructure implementing the ring of integers modulo a prime p
*/

template <typename int p>
class ModP {
    private:
        int value = 0;

        inline int mod(int x) const {
            while (x < 0) x += p;
			x %= p;
			assert(0 <= x && x < p);
            return x;
        }

    public:

	//-------------------------------------------------------------------------------------------------------------|
	// Constructors
	//-------------------------------------------------------------------------------------------------------------|
	ModP() = default;
	ModP(int val) : value(mod(val)) {};

	//-------------------------------------------------------------------------------------------------------------|
	// Operators
	//-------------------------------------------------------------------------------------------------------------|
	ModP operator+=(const ModP& rhs) { value = mod(value + rhs.value); return *this; }
	ModP operator-=(const ModP& rhs) { value = mod(value - rhs.value); return *this; }
	ModP operator*=(const ModP& rhs) { value = mod(value * rhs.value); return *this; }
	ModP operator%=(const ModP& rhs) { value = 0; 					   return *this; }
	ModP operator/=(const ModP& rhs) { // Brute force, kinda ugly
        if (rhs.value == 0) {
			throw std::invalid_argument("Cannot divide by zero.");
		}
        for (int i = 0; i < p; ++i) {
            if (i * rhs.value == value) {
                value = i;
                return *this;
            }
        }
		// Will only occur iff p is not prime and rhs a divisor of p
        throw std::invalid_argument("Not a divisor.");
	}

	friend ModP operator+(ModP lhs, const ModP& rhs) { return lhs += rhs; }
	friend ModP operator-(ModP lhs, const ModP& rhs) { return lhs -= rhs; }
	friend ModP operator*(ModP lhs, const ModP& rhs) { return lhs *= rhs; }
	friend ModP operator/(ModP lhs, const ModP& rhs) { return lhs /= rhs; }
	friend ModP operator%(ModP lhs, const ModP& rhs) { return lhs %= rhs; }

	bool operator==(const ModP& rhs) const { return value == rhs.value; }
	bool operator!=(const ModP& rhs) const { return !(*this == rhs); }

	//-------------------------------------------------------------------------------------------------------------|
	// Getters
	//-------------------------------------------------------------------------------------------------------------|
	int getValue() const { return value; }
};

template<typename int p>
std::ostream& operator<<(std::ostream& os, const ModP<p>& obj) { 
	return os << "[" << obj.getValue() << "]";
}