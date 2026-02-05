#pragma once
#include <cassert>
#include <ostream>
#include <sstream>
#include <vector>

#include "util/helper.hpp"

template <typename R> class Polynomial;

template <typename R> std::ostream& operator<<(std::ostream&, const Polynomial<R>&);
template <typename R> void swap(Polynomial<R>&, Polynomial<R>&);

/**
* @brief A class representing a polynomial with coefficients in a commutative ring R.
*
* @tparam R The type of the underlying commutative ring
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

	R* coeffs = nullptr;	///< Pointer to the array containing the (degree + 1) coefficients
	int degree = -1;		///< The degree

	//---------------------------------------------------------------------|
	// Polynomial division
	//---------------------------------------------------------------------|

	bool polyDiv(const Polynomial&, Polynomial&, Polynomial&) const;

	//---------------------------------------------------------------------|
	// Memory access
	//---------------------------------------------------------------------|

	Polynomial(R* arr, int deg) : degree(deg), coeffs(arr) {};	///< Standard constructor

//-------------------------------------------------------------------------------------------------------------|
// Public
//-------------------------------------------------------------------------------------------------------------|
public:
	//---------------------------------------------------------------------|
	// Constructors
	//---------------------------------------------------------------------|

	Polynomial() = default;
	Polynomial(const R&, int = 0);				///< Constant polynomial constructor
	Polynomial(int);							///< Monomial constructor
	Polynomial(const Polynomial&);				///< Copy constructor
	Polynomial(Polynomial&&) noexcept;			///< Move constructor

	~Polynomial() { util::deallocate(coeffs); }	///< Destructor

	//---------------------------------------------------------------------|
	// Operators (in place)
	//---------------------------------------------------------------------|

	Polynomial operator+=(const Polynomial& rhs) { *this = *this + rhs; return *this; }	///< In place polynomial addition
	Polynomial operator-=(const Polynomial& rhs) { *this = *this - rhs; return *this; }	///< In place polynomial subtraction
	Polynomial operator*=(const Polynomial& rhs) { *this = *this * rhs; return *this; }	///< In place polynomial multiplication
	Polynomial operator/=(const Polynomial& rhs) { *this = *this / rhs; return *this; } ///< In place polynomial division
	Polynomial operator%=(const Polynomial& rhs) { *this = *this % rhs; return *this; }	///< In place polynomial modulation

	Polynomial& operator=(const Polynomial&);											///< Assignment operator
	Polynomial& operator=(Polynomial&&) noexcept;										///< Move operator

	//---------------------------------------------------------------------|
	// Operators (out of place)
	//---------------------------------------------------------------------|

	Polynomial operator+(const Polynomial&) const;								///< Out of place polynomial addition
	Polynomial operator-(const Polynomial&) const;								///< Out of place polynomial subtraction
	Polynomial operator*(const Polynomial&) const;								///< Out of place polynomial multiplication
	Polynomial operator/(const Polynomial&) const;
	Polynomial operator%(const Polynomial&) const;
	
	bool operator==(const Polynomial&) const;									///< Polynomial equality
	bool operator!=(const Polynomial& rhs) const { return !(*this == rhs); }	///< Polynomial inequality

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

	inline int getDegree() const { return degree; };	///< Getter for the degree

	//---------------------------------------------------------------------|
	// Others
	//---------------------------------------------------------------------|

	friend void swap<>(Polynomial&, Polynomial&);		///< Swaps two polynomials

	//---------------------------------------------------------------------|
	// Static functions
	//---------------------------------------------------------------------|

	static bool parseToPolynomial(const std::string&, Polynomial&);
	static int characteristic() { return R::characteristic(); }		///< Returns the characteristic of R[X]
};

#include "polynomial.ipp"