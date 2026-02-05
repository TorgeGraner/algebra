#pragma once
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

#include "util/helper.hpp"

#include "datastructure/polynomial.hpp"

template <typename R> class Matrix;
template <typename R> void swap(Matrix<R>&, Matrix<R>&);

/**
 * @brief A class representing a matrix with entries in a ring R
 * 
 * @tparam R The type of the entries of the matrix (for some functions required to be a euclidean ring)
 * 
 * @todo 1) Implement L-R Decomposition
 * @todo 2) Representation for zero-blocks or diagonal/triangular/block/scalar matrices
 * @todo 3) Minimal polynomial for frobenius normal form
 */
template<typename R>
class Matrix {
template <typename T> friend class Matrix;
//-------------------------------------------------------------------------------------------------------------|
// Private
//-------------------------------------------------------------------------------------------------------------|
private:
	//---------------------------------------------------------------------|
	// Member variables
	//---------------------------------------------------------------------|

	R* entries = nullptr;	///< Pointer to the array containing the entries

	int n = -1; 			///< Number of rows
	int m = -1; 			///< Number of columns
	
	bool reduced = false;	///< True if matrix is in rref
	
	R determinant;			///< The determinant of the matrix
	int rank = -1;			///< The rank of the matrix

	//---------------------------------------------------------------------|
	// Memory access
	//---------------------------------------------------------------------|

	Matrix(R* arr, int n, int m) : entries(arr), n(n), m(m), reduced(false) {}	///< Standard constructor

	inline int index(int i, int j, int m) const { return m * i + j; }			///< Converts tuple to array index
	inline R& rawEntry(int i, int j) const { return entries[index(i, j, m)]; }	///< Returns reference to entry (i, j)
	
	//---------------------------------------------------------------------|
	// Elementary row operations (in place)
	//---------------------------------------------------------------------|

	void swapRows(int, int);				///< Swaps two rows
	void scaleRow(int, const R&);			///< Scales a row by a constant
	void cancelRow(int, const R&);			///< Cancels a divisor from a row
	void subtractRow(int, int, const R&);	///< Subtracts a multiple of one row from another

	//---------------------------------------------------------------------|
	// Row reduction
	//---------------------------------------------------------------------|

	void rref(bool);

	//---------------------------------------------------------------------|
	// Characteristic polynomial algorithms
	//---------------------------------------------------------------------|

	Polynomial<R> charPolyNaive() const;
	Polynomial<R> faddeevLeVerrier() const;
	
//-------------------------------------------------------------------------------------------------------------|
// Public
//-------------------------------------------------------------------------------------------------------------|
public:
	//------------------------------------------------------------------------------|
	// Constructors
	//------------------------------------------------------------------------------|

	Matrix() = default;							///< Default constructor
	Matrix(R, const int, const int);			///< Constant matrix constructor
	Matrix(const Matrix&);						///< Copy constructor
	Matrix(Matrix&&) noexcept;					///< Move constructor
	
	~Matrix() { util::deallocate(entries); }	///< Destructor

	//------------------------------------------------------------------------------|
	// Operators (in place)
	//------------------------------------------------------------------------------|

	Matrix operator+=(const Matrix&);												///< In place matrix addition
	Matrix operator-=(const Matrix&);												///< In place matrix subtraction						
	Matrix operator*=(const Matrix& rhs) { *this = *this * rhs; return *this; }		///< In place matrix multiplication
	Matrix operator*=(const R&);													///< In place scalar multiplication
	
	Matrix& operator=(const Matrix&);												///< Assignment operator
	Matrix& operator=(Matrix&&) noexcept;											///< Move operator

	//-------------------------------------------------------------------------------|
	// Operators (out of place)
	//-------------------------------------------------------------------------------|

	friend Matrix operator+(Matrix lhs, const Matrix& rhs) { return lhs += rhs; }	///< Out of place matrix addition
	friend Matrix operator-(Matrix lhs, const Matrix& rhs) { return lhs -= rhs; }	///< Out of place matrix subtraction
	friend Matrix operator*(Matrix lhs, const R& rhs)	    { return lhs *= rhs; }	///< Out of place scalar multiplication
	Matrix operator*(const Matrix&) const;											///< Out of place matrix multiplication

	bool operator==(const Matrix<R>&) 	  const;									///< Matrix equality
	bool operator==(const R&) 			  const;									///< Matrix equality
	bool operator!=(const Matrix<R>& rhs) const { return !(*this == rhs); }			///< Matrix inequality
	bool operator!=(const R& rhs)         const { return !(*this == rhs); }			///< Scalar inequality

	const R& operator()(const int, const int) const;
	
	friend std::ostream& operator<< <>(std::ostream&, const Matrix&);

	//-------------------------------------------------------------------------------------------------------------|
	// Operations
	//-------------------------------------------------------------------------------------------------------------|

	R inverse(Matrix&) const;

	Matrix nullspace() const;
	Matrix complete(const Matrix&) const;
	Matrix mergeHorizontal(const Matrix&) const;

	Polynomial<R> charPoly() const;

	//-------------------------------------------------------------------------------------------------------------|
	// Getters
	//-------------------------------------------------------------------------------------------------------------|

	inline int getN() const { return n; }	///< Getter for the number of rows
	inline int getM() const { return m; }	///< Getter for the number of columns

	int getRank() const;					///< Getter for the rank
	R getDeterminant() const;				///< Getter for the determinant
	
	Matrix getRows(int, int) const;
	Matrix getColumns(int, int) const;

	//-------------------------------------------------------------------------------------------------------------|
	// Others
	//-------------------------------------------------------------------------------------------------------------|

	friend void swap<>(Matrix&, Matrix&);	///< Swaps two matrices

	//-------------------------------------------------------------------------------------------------------------|
	// Static functions
	//-------------------------------------------------------------------------------------------------------------|

	static bool parseToMatrix(const std::string&, Matrix&, int = -1);

	/**
	 * @brief Converts this matrix to a matrix of another type
	 * 
	 * @tparam R_OUT The target ring type
	 * 
	 * @return The converted matrix
	 */
	template<typename R_OUT>
	Matrix<R_OUT> convert() const { return Matrix<R_OUT>(util::copy<R, R_OUT>(entries, n * m), n, m); }

	static int characteristic() { return R::characteristic(); }		///< Returns the characteristic of @f[R^{n\times m}@f]
};

#include "matrix.ipp"