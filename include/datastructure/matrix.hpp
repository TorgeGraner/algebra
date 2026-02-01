#pragma once
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "util/helper.hpp"

#include "datastructure/polynomial.hpp"

template <typename R> class Matrix;
template <typename R> void swap(Matrix<R>&, Matrix<R>&);

/* 
* A class representing a matrix over an integral domain or a field R, no multiplicative inverses are generally required,
* however if R is not a field, the results of some functions as the inverse matrix will return scaled versions
*
* Note: Always include this file AFTER the inclusion of the class implementing R
* 
* TODO: 
* 1) Implement L-R Decomposition
* 2) Representation for zero-blocks or diagonal/triangular/block/scalar matrices (maybe some sort of quads?)
* 3) Minimal polynomial for frobenius normal form
*/
template<typename R>
class Matrix {
template <typename T> friend class Matrix; // Required for convert function
//-------------------------------------------------------------------------------------------------------------|
// Private
//-------------------------------------------------------------------------------------------------------------|
private:
	//---------------------------------------------------------------------|
	// Member variables
	//---------------------------------------------------------------------|

	R* entries = nullptr;	// Pointer to the array containing the entries

	int n = -1; 			// Number of rows
	int m = -1; 			// Number of columns
	
	bool reduced = false;	// True if matrix is in rref
	
	R determinant;			// The determinant of the matrix
	int rank = -1;			// The rank of the matrix

	//---------------------------------------------------------------------|
	// Memory access
	//---------------------------------------------------------------------|

	Matrix(R* arr, int n, int m) : entries(arr), n(n), m(m), reduced(false) {}	// Standard constructor

	inline int index(int i, int j, int m) const { return m * i + j; }			// Converts tuple to array index
	inline R& rawEntry(int i, int j) const { return entries[index(i, j, m)]; };	// Returns reference to entry (i, j)
	
	//---------------------------------------------------------------------|
	// Elementary row operations (in place)
	//---------------------------------------------------------------------|

	void swapRows(int, int);
	void scaleRow(int, const R&);
	void cancelRow(int, const R&);
	void subtractRow(int, int, const R&);

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

	Matrix() = default;							// Default constructor
	Matrix(R, const int, const int);
	Matrix(const Matrix&);
	Matrix(Matrix&&) noexcept;
	
	~Matrix() { util::deallocate(entries); }	// Destructor

	//------------------------------------------------------------------------------|
	// Operators (in place)
	//------------------------------------------------------------------------------|

	Matrix operator+=(const Matrix&);
	Matrix operator-=(const Matrix&);
	Matrix operator*=(const Matrix& rhs) { *this = *this * rhs; return *this; }
	Matrix operator*=(const R&);
	
	Matrix& operator=(const Matrix&);
	Matrix& operator=(Matrix&&) noexcept;

	//-------------------------------------------------------------------------------|
	// Operators (out of place)
	//-------------------------------------------------------------------------------|

	friend Matrix operator+(Matrix lhs, const Matrix& rhs) { return lhs += rhs; }
	friend Matrix operator-(Matrix lhs, const Matrix& rhs) { return lhs -= rhs; }
	friend Matrix operator*(Matrix lhs, const R& rhs)	   { return lhs *= rhs; }
	Matrix operator*(const Matrix&) const;

	bool operator==(const Matrix<R>&) 	  const;
	bool operator==(const R&) 			  const;
	bool operator!=(const Matrix<R>& rhs) const { return !(*this == rhs); }
	bool operator!=(const R& rhs)         const { return !(*this == rhs); }

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

	inline int getN() const { return n; }
	inline int getM() const { return m; }

	int getRank() const;
	R getDeterminant() const;
	
	Matrix getRows(int, int) const;
	Matrix getColumns(int, int) const;

	//-------------------------------------------------------------------------------------------------------------|
	// Others
	//-------------------------------------------------------------------------------------------------------------|

	friend void swap<>(Matrix&, Matrix&);

	//-------------------------------------------------------------------------------------------------------------|
	// Static functions
	//-------------------------------------------------------------------------------------------------------------|

	static bool parseToMatrix(const std::string&, Matrix&, int = -1);

	// Converts matrix to another ring
	template<typename R_OUT>
	Matrix<R_OUT> convert() const {
		return Matrix<R_OUT>(util::copy<R, R_OUT>(entries, n * m), n, m);
	}
};

#include "matrix.tpp"