#pragma once
#include <ostream>

template <typename int epsilon> class EpsFloat;
template <typename int epsilon> std::ostream& operator<<(std::ostream&, const EpsFloat<epsilon>&);

// A datastructure implementing the ring of integers modulo p, if p is prime this is a field
template <typename int epsilon>
class EpsFloat {
    private:
        double value = 0;

    public:

	//-------------------------------------------------------------------------------------------------------------|
	// Constructors
	//-------------------------------------------------------------------------------------------------------------|

	EpsFloat() = default;
	EpsFloat(double val) : value(val) {};

	//-------------------------------------------------------------------------------------------------------------|
	// Operators
	//-------------------------------------------------------------------------------------------------------------|

	EpsFloat operator+=(const EpsFloat& rhs) { value += rhs.value; return *this; }
	EpsFloat operator-=(const EpsFloat& rhs) { value -= rhs.value; return *this; }
	EpsFloat operator*=(const EpsFloat& rhs) { value *= rhs.value; return *this; }
	EpsFloat operator/=(const EpsFloat& rhs) { value /= rhs.value; return *this; }
	EpsFloat operator%=(const EpsFloat& rhs) { value = 0;          return *this; }

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