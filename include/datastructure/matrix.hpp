#pragma once
#include <iostream>
#include <cassert>

#include "util/helper.hpp"

template <typename R> class Matrix;
template <typename R> void swap(Matrix<R>&, Matrix<R>&);

/* 
* @brief A class representing a matrix over a given Ring R, no multiplicative inverses are generally required,
* however if R is not a field, the results of some functions as the inverse matrix will return scaled versions
* 
* TODO: 
* 1) Implement L-R Decomposition, change calcRef
* 2) Representation for zero-blocks or diagonal/triangular/block/scalar matrices (maybe some sort of quads?)
* 3) Minimal polynomial for frobenius normal form
* 4) Improve Print operator
*/
template<typename R>
class Matrix {
private:
	R* entries = nullptr;
	int n = -1; 
	int m = -1;
	
public:
	//-------------------------------------------------------------------------------------------------------------|
	// Constructors
	//-------------------------------------------------------------------------------------------------------------|
	Matrix() = default;							// Default constructor
	Matrix(R*, const int, const int);			// Standard constructor
	Matrix(R, const int, const int m);			// Diagonal matrix constructor
	Matrix(const Matrix&);						// Copy constructor
	Matrix(Matrix&&) noexcept;					// Move constructor

	~Matrix() { util::deallocate(entries); };	// Destructor
	//-------------------------------------------------------------------------------------------------------------|
	// Operators
	//-------------------------------------------------------------------------------------------------------------|
	Matrix operator+=(const Matrix&);
	Matrix operator-=(const Matrix&);
	Matrix operator*=(const Matrix&);

	friend Matrix operator+(Matrix lhs, const Matrix& rhs) { return lhs += rhs; }
	friend Matrix operator-(Matrix lhs, const Matrix& rhs) { return lhs -= rhs; }
	friend Matrix operator*(Matrix lhs, const Matrix& rhs) { return lhs *= rhs; }

	bool operator==(const Matrix<R>&) const;
	bool operator==(const R&) const;
	bool operator!=(const Matrix<R>& rhs) const { return !(*this == rhs); }
	bool operator!=(const R& rhs) const { return !(*this == rhs); }

	R& operator()(const int, const int) const;

	// Functions required since the destructor is implemented
	Matrix& operator=(const Matrix&);		// Assignment operator
	Matrix& operator=(Matrix&&) noexcept;	// Move operator

	friend void swap<>(Matrix&, Matrix&);

	friend std::ostream& operator<< <>(std::ostream&, const Matrix&);

	//-------------------------------------------------------------------------------------------------------------|
	// Getters
	//-------------------------------------------------------------------------------------------------------------|
	const int getN() const { return n; }
	const int getM() const { return m; }
};

//-------------------------------------------------------------------------------------------------------------|
// Constructors
//-------------------------------------------------------------------------------------------------------------|
// Standard constructor
template <typename R>
Matrix<R>::Matrix(R* _entries, const int n, const int m) : n(n), m(m) {
	if (n * m != 0) {
		entries = util::allocate<R>(n * m);
		std::copy(_entries, _entries + n * m, entries);
	}
}

// Diagonal quadratic matrix constructor
template <typename R>
Matrix<R>::Matrix(R entry, const int n, int m) : n(n), m(m) {
	entries = util::allocate<R>(n * m);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			if (i == j) {
				entries[m * i + j] = entry;
			}
			else {
				entries[m * i + j] = 0;
			}
		}
	}
}

// Copy Constructor
template<typename R>
Matrix<R>::Matrix(const Matrix<R>& orig) : n(orig.n), m(orig.m) {
	entries = util::allocate<R>(n * m);
	std::copy(orig.entries, orig.entries + n * m, entries);
}

// Move constructor
template <typename R>
Matrix<R>::Matrix(Matrix<R>&& src) noexcept : Matrix{} {
	swap(*this, src);
}

//-------------------------------------------------------------------------------------------------------------|
// Operators
//-------------------------------------------------------------------------------------------------------------|
template<typename R>
Matrix<R> Matrix<R>::operator+=(const Matrix<R>& rhs) {
	if (m != rhs.m && n != rhs.n) {
		throw std::invalid_argument("Cannot add matrices of different dimensions.");
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			entries[m * i + j] += rhs(i, j);
		}
	}
	return *this;
}

template<typename R>
Matrix<R> Matrix<R>::operator-=(const Matrix<R>& rhs) {
	if (m != rhs.m && n != rhs.n) {
		throw std::invalid_argument("Cannot subtract matrices of different dimensions.");
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			entries[m * i + j] -= rhs(i, j);
		}
	}
	return *this;
}

// Matrix multiplication
template<typename R>
Matrix<R> Matrix<R>::operator*=(const Matrix<R>& multiplicand) {
	int multM = multiplicand.getM();
	if (multiplicand.getN() == 1 && multM == 1) {
		for (int k = 0; k < n * m; ++k) {
			entries[k] *= multiplicand(0, 0);
		}
		return *this;
	}
	if (m != multiplicand.getN()) {
		throw std::invalid_argument("Cannot multiply matrices with wrong dimensions.");
	}
	R* arr = util::allocate<R>(n * multM);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < multM; ++j) {
			int h = i * multM + j;
			arr[h] = 0;
			for (int k = 0; k < m; ++k) {
				arr[h] += (*this)(i, k) * multiplicand(k, j);
			}
		}
	}
	m = multM;
	util::deallocate<R>(entries);
	entries = arr;
	return *this;
}

// Comparison to another matrix
template <typename R>
bool Matrix<R>::operator==(const Matrix<R>& rhs) const {
	if (n != rhs.n || m != rhs.m) {
		return false;
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			if ((*this)(i, j) != rhs(i, j)) {
				return false;
			}
		}
	}
	return true;
}

// Check if matrix is scaled identity matrix
template <typename R>
bool Matrix<R>::operator==(const R& rhs) const {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			if (i == j) {
				if ((*this)(i, j) != rhs) return false;
			}
			else {
				if ((*this)(i, j) != 0) return false;
			}
		}
	}
	return true;
}

// Assignment operator
template <typename R>
Matrix<R>& Matrix<R>::operator=(const Matrix<R>& rhs) {
	Matrix tmp(rhs);
	swap(*this, tmp);
	return *this;
}

// Move operator
template <typename R>
Matrix<R>& Matrix<R>::operator=(Matrix<R>&& src) noexcept {
	Matrix tmp(src);
	swap(*this, src);
	return *this;
}

// Getter for j-th entry in i-th row
template <typename R>
R& Matrix<R>::operator()(const int i, const int j) const {
	if (0 > i || i >= n || 0 > j || j >= m) {
		throw std::invalid_argument("Index out of bounds for matrix.");
	}
	return entries[m * i + j];
}

// Printing operator
template <typename R>
std::ostream& operator<< <>(std::ostream& os, const Matrix<R>& obj) {
	std::vector<int> colSizes;
	std::ostringstream helpOs; 
	// Calculate max size of spaces in each column
	for (int j = 0; j < obj.m; ++j) {
		int mSize = 0;
		for (int i = 0; i < obj.n; ++i) {
			auto start = helpOs.tellp();
			helpOs << obj(i, j);
			auto end = helpOs.tellp();
			mSize = std::max(mSize, int(end - start));
			helpOs.clear();
		}
		colSizes.push_back(mSize);
	}
	// Print and fill
	for (int i = 0; i < obj.n; ++i) {
		for (int j = 0; j < obj.m; ++j) {
			auto start = helpOs.tellp();
			helpOs << obj(i, j);
			auto end = helpOs.tellp();

			os << obj(i, j);
			int fill = colSizes[j] - int(end - start) + 1;
			for (int k = 0; k < fill; ++k) os << " ";
		}
		os.flush();
		os << "\n";
	}
	return os;
}

template <typename R>
void swap(Matrix<R>& lhs, Matrix<R>& rhs) {
	std::swap(lhs.n, rhs.n);
	std::swap(lhs.m, rhs.m);
	std::swap(lhs.entries, rhs.entries);
}
