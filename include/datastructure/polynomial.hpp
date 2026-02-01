#pragma once
#include <cassert>
#include <ostream>
#include <sstream>
#include <vector>

#include "util/helper.hpp"

template <typename R> class Polynomial;

template <typename R> std::ostream& operator<<(std::ostream&, const Polynomial<R>&);
template <typename R> void swap(Polynomial<R>&, Polynomial<R>&);

/* 
* A class representing a polynomial with coefficients in a ring R.
*
* Note: Always include this file AFTER the inclusion of the class implementing R
*/
template <typename R> 
class Polynomial {
template<typename U> friend class Matrix;
//-------------------------------------------------------------------------------------------------------------|
// Private
//-------------------------------------------------------------------------------------------------------------|
private:
	//---------------------------------------------------------------------|
	// Member variables
	//---------------------------------------------------------------------|

	R* coeffs = nullptr;	// Pointer to the array containing the (degree + 1) coefficients
	int degree = -1;		// The degree

	//---------------------------------------------------------------------|
	// Polynomial division
	//---------------------------------------------------------------------|

	bool polyDiv(const Polynomial&, Polynomial&, Polynomial&) const;

	//---------------------------------------------------------------------|
	// Memory access
	//---------------------------------------------------------------------|

	Polynomial(R* arr, int deg) : degree(deg), coeffs(arr) {};	// Standard constructor

//-------------------------------------------------------------------------------------------------------------|
// Public
//-------------------------------------------------------------------------------------------------------------|
public:
	//---------------------------------------------------------------------|
	// Constructors
	//---------------------------------------------------------------------|

	Polynomial() = default;						// Default constructor
	Polynomial(const R&, int = 0);
	Polynomial(int);
	Polynomial(const Polynomial&);
	Polynomial(Polynomial&&) noexcept;

	~Polynomial() { util::deallocate(coeffs); }	// Destructor

	//---------------------------------------------------------------------|
	// Operators (in place)
	//---------------------------------------------------------------------|

	Polynomial operator+=(const Polynomial& rhs) { *this = *this + rhs; return *this; }
	Polynomial operator-=(const Polynomial& rhs) { *this = *this - rhs; return *this; }
	Polynomial operator*=(const Polynomial& rhs) { *this = *this * rhs; return *this; }
	Polynomial operator/=(const Polynomial& rhs) { *this = *this / rhs; return *this; }
	Polynomial operator%=(const Polynomial& rhs) { *this = *this % rhs; return *this; }

	Polynomial& operator=(const Polynomial&);
	Polynomial& operator=(Polynomial&&) noexcept;

	//---------------------------------------------------------------------|
	// Operators (out of place)
	//---------------------------------------------------------------------|

	Polynomial operator+(const Polynomial&) const;
	Polynomial operator-(const Polynomial&) const;
	Polynomial operator*(const Polynomial&) const;
	Polynomial operator/(const Polynomial&) const;
	Polynomial operator%(const Polynomial&) const;
	
	bool operator==(const Polynomial&) const;
	bool operator!=(const Polynomial& rhs) const { return !(*this == rhs); }

	const R& operator()(int) const;
	
	friend std::ostream& operator<< <>(std::ostream&, const Polynomial&);

	//---------------------------------------------------------------------|
	// Operations
	//---------------------------------------------------------------------|

	R map(const R&) const;
	Polynomial derivative() const;

	//---------------------------------------------------------------------|
	// Getters
	//---------------------------------------------------------------------|

	inline int getDegree() const { return degree; };

	//---------------------------------------------------------------------|
	// Others
	//---------------------------------------------------------------------|

	friend void swap<>(Polynomial&, Polynomial&);

	//---------------------------------------------------------------------|
	// Static functions
	//---------------------------------------------------------------------|

	static bool parseToPolynomial(const std::string&, Polynomial&);
	static int characteristic() { return R::characteristic(); }		// Returns the characteristic of R[X]
};

#include "polynomial.tpp"