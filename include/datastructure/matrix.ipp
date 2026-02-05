#pragma once
//-------------------------------------------------------------------------------------------------------------|
// PRIVATE
//---------------------------------------------------------------------|---------------------------------------|
// Elementary row operations (in place)
//---------------------------------------------------------------------|

template<typename R>
void Matrix<R>::scaleRow(int i, const R& scale) {
	for (int j = 0; j < m; ++j) rawEntry(i, j) *= scale;
}

template<typename R>
void Matrix<R>::cancelRow(int i, const R& div) {
	for (int j = 0; j < m; ++j) rawEntry(i, j) /= div;
}

template<typename R> 
void Matrix<R>::swapRows(int i1, int i2) {
	for (int j = 0; j < m; ++j) std::swap(rawEntry(i1, j), rawEntry(i2, j));
}

template<typename R>
void Matrix<R>::subtractRow(int i1, int i2, const R& scale) {
	for (int j = 0; j < m; ++j) rawEntry(i1, j) -= scale * (*this)(i2, j);
}

//---------------------------------------------------------------------|
// Row reduction
//---------------------------------------------------------------------|

/**
 * @brief Turns the matrix into (reduced) row echelon form using the Gaussian algorithm
 *
 * @tparam R The type of the underlying euclidean ring
 * @param reduced Indicates if the row echelon form should be reduced (rref, with zeros above the diagonal) or not (ref)
 * 
 * @note Sets the member variables rank, reduced and determinant, if applicable
 */
template <typename R>
void Matrix<R>::rref(bool reduced) {
	R detNumer = 1;
	R detDenom = 1;
	rank = 0;
	// The current number of reduced columns
	int numReduced = 0;
	// The row index of the current non-zero element
	int nzRowInd = -1;

	for (int j = 0; j < m; ++j) {
		nzRowInd = -1;
		// Find the first non-zero entry and its row-index s in column j
		for (int i = numReduced; i < n; ++i) {
			if ((*this)(i, j) != 0) {
				nzRowInd = i;
				break;
			}
		}
		// If there is no non zero entry, move to the next column
		if (nzRowInd == -1) {
			detNumer = 0;
			continue;
		}
		++rank;
		// Eradicate all other non-zero entries in column j
		int lower = reduced ? 0 : numReduced;
		for (int i = lower; i < n; ++i) {
			if (i != nzRowInd && (*this)(i, j) != 0) {
				R gcd = util::gcd((*this)(nzRowInd, j), (*this)(i, j));

				R nzRowFact = (*this)(i, j) / gcd;
				R elRowFact = (*this)(nzRowInd, j) / gcd;
				// Reduce row i and keep it as simple as possible
				scaleRow(i, elRowFact);
				subtractRow(i, nzRowInd, nzRowFact);
				detDenom *= elRowFact;
			}
			// Calculate the gcd of row i and factor it out
			R rowGcd = 0;
			for (int j = 0; j < m; ++j) {
				rowGcd = util::gcd(rowGcd, (*this)(i, j));
			}
			if (rowGcd != 0 && rowGcd != 1) {
				cancelRow(i, rowGcd);
				detNumer *= rowGcd;
			}
		}
		// Switch rows nzRowInd and num
		if (nzRowInd != numReduced) {
			swapRows(nzRowInd, numReduced);
			detNumer *= (-1);
		}
		++numReduced;
	}
	if (n == m) {
		for (int i = 0; i < n; ++i) {
			detNumer *= (*this)(i, i);
		}
		determinant = detNumer / detDenom;
	}
	this->reduced = true;
}

//---------------------------------------------------------------------|
// Characteristic polynomial algorithms
//---------------------------------------------------------------------|

/**
 * @brief Calculate the characteristic polynomial using the naive definition @f[ p(X) := det(X*I - A) @f]
 * 
 * @tparam R The type of the underlying euclidean ring
 * 
 * @warning Highly volatile and should only be used with rings of finite characteristic
 * 
 * @return The characteristic polynomial of the matrix
 */
template <typename R>
Polynomial<R> Matrix<R>::charPolyNaive() const {
	Matrix<Polynomial<R>>polyMat(Polynomial<R>(1, 1), n, m);
	polyMat -= convert<Polynomial<R>>();
	return polyMat.getDeterminant();
}

/**
 * @brief Calculate the characteristic polynomial using the faddeevLeVerrier algorithm
 * 
 * @tparam R The type of the underlying ring
 * 
 * @warning Should only be used with rings of characteristic zero
 * 
 * @return The characteristic polynomial of the matrix
 */
template<typename R>
Polynomial<R> Matrix<R>::faddeevLeVerrier() const {
	Matrix<R> C(*this);
	R* arr = util::allocate<R>(n + 1);
	arr[n] = 1;
	for (int k = 1; k <= n; ++k) {
		if (k > 1) {
			C += Matrix<R>(arr[n - k + 1], n, n);
			C = (*this) * C; // Note how *= is right sided multiplication!
		}
		arr[n - k] = 0;
		for (int j = 0; j < n; ++j) {
			arr[n - k] += C(j, j);
		}
		// Possible error in finite rings
		arr[n - k] /= (-k);
	}
	return Polynomial<R>(arr, n);
}

//-------------------------------------------------------------------------------------------------------------|
// PUBLIC
//----------------------------------------------------------------------|--------------------------------------|
// Constructors
//----------------------------------------------------------------------|

template <typename R>
Matrix<R>::Matrix(R entry, const int n, int m) : n(n), m(m), reduced(false) {
	entries = util::zeroes<R>(n * m);
	for (int i = 0; i < std::min(n, m); ++i) {
		rawEntry(i, i) = entry;
	}
}

template<typename R>
Matrix<R>::Matrix(const Matrix<R>& orig) : n(orig.n), m(orig.m), reduced(orig.reduced) {
	entries = util::copy<R, R>(orig.entries, n * m);
}

template <typename R>
Matrix<R>::Matrix(Matrix<R>&& src) noexcept : Matrix{} {
	swap(*this, src);
}

//----------------------------------------------------------------------|
// Operators (in place)
//----------------------------------------------------------------------|

template<typename R>
Matrix<R> Matrix<R>::operator+=(const Matrix<R>& rhs) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			rawEntry(i, j) += rhs(i, j);
		}
	}
	return *this;
}

template<typename R>
Matrix<R> Matrix<R>::operator-=(const Matrix<R>& rhs) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			rawEntry(i, j) -= rhs(i, j);
		}
	}
	return *this;
}

template<typename R>
Matrix<R> Matrix<R>::operator*=(const R& rhs) {
	for (int k = 0; k < n * m; ++k) entries[k] *= rhs;
	return *this;
}

template <typename R>
Matrix<R>& Matrix<R>::operator=(const Matrix<R>& rhs) {
	Matrix tmp(rhs);
	swap(*this, tmp);
	return *this;
}

template <typename R>
Matrix<R>& Matrix<R>::operator=(Matrix<R>&& src) noexcept {
	Matrix tmp(src);
	swap(*this, src);
	return *this;
}

//-------------------------------------------------------------------------------|
// Operators (out of place)
//-------------------------------------------------------------------------------|

template<typename R>
Matrix<R> Matrix<R>::operator*(const Matrix<R>& multiplicand) const {
	int multM = multiplicand.getM();
	R* arr = util::allocate<R>(n * multM);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < multM; ++j) {
			int id = index(i, j, multM);
			arr[id] = 0;
			for (int k = 0; k < m; ++k) {
				arr[id] += (*this)(i, k) * multiplicand(k, j);
			}
		}
	}
	return Matrix<R>(arr, n, multM);
}

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

template <typename R>
bool Matrix<R>::operator==(const R& rhs) const {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			if (i == j) {
				if ((*this)(i, j) != rhs) return false;
			} else {
				if ((*this)(i, j) != 0) return false;
			}
		}
	}
	return true;
}

/**
 * @brief Getter for entry (i, j)
 * 
 * @tparam R The underlying type
 * @param i The row index
 * @param j The column index
 * 
 * @throws std::invalid_argument if the indices are out of bounds
 * 
 * @return The entry (i, j)
 */
template <typename R>
const R& Matrix<R>::operator()(const int i, const int j) const {
	if (0 > i || i >= n || 0 > j || j >= m) {
		throw std::invalid_argument("Index out of bounds for matrix.");
	}
	return entries[index(i, j, m)];
}

/**
 * @brief Stream operator
 * 
 * @tparam R The underlying type
 * @param os The output stream
 * @param obj The object to stream
 * @return The output stream
 */
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

//----------------------------------------------------------------------|
// Operations
//----------------------------------------------------------------------|

/**
 * @brief Calculate the minimally scaled inverse of the matrix 
 * 
 * @tparam R The type of the underlying euclidean ring
 * 
 * @param[out] inv The inverse matrix scaled by a factor
 * 
 * @return R The scaling factor such that @f[ inv = scale * A^{-1} @f]
 */
template <typename R>
R Matrix<R>::inverse(Matrix<R>& inv) const {
	// Get (mat, I)
	Matrix<R> reduced = mergeHorizontal(Matrix<R>(1, n, n));
	// Reduce, resulting matrix contains inverse multiplied by a diagonal matrix
	reduced.rref(true);
	// Find least common multiple of elements 
	R diagLcm = 1;
	for (int i = 0; i < n; ++i) {
		R temp = reduced(i, i);
		if (temp == 0) {
			throw std::invalid_argument("Matrix is not invertible.");
		}
		diagLcm *= temp / util::gcd(diagLcm, temp);
	}
	// Scale mat^(-1) according to the lcm and the entry in the diagonal matrix in the corresponding row
	R* arr = util::allocate<R>(n * m);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			arr[i * m + j] = reduced(i, m + j) * (diagLcm / reduced(i, i));
		}
	}
	inv = Matrix<R>(arr, n, m);
	return diagLcm;
}

/**
 * @brief Return a base for the nullspace of this matrix
 * 
 * @tparam R The type of the underlying euclidean ring
 * 
 * @return A matrix whose columns form a bas of the nullspace of this matrix
 */
template<typename R>
Matrix<R> Matrix<R>::nullspace() const {
	Matrix<R> phi(*this);
	phi.rref(true);
	// Row indicates the current row in the rref, column counts the current number of vectors in the incomplete base
	int row = 0;
	int column = 0;
	int dim = m - phi.getRank();
	// The indices of the first non-zero element in the rref
	std::vector<int> indices(phi.getRank());
	R* arr = util::zeroes<R>(m * dim);

	for (int j = 0; j < m; ++j) {
		// if A[row][j] == 0, solve for variable x_j
		if (row >= n || phi(row, j) == 0) {
			// Calculate LCM of non-zero elements in column to avoid multiplicative inverses and unnecessarily big entries
			R colLcm = 1;
			for (int k = 0; k < row; ++k) {
				if (phi(k, j) != 0) {
					colLcm *= phi(k, indices[k]) / util::gcd(colLcm, phi(k, indices[k]));
				}
			}
			// Add the vector corresponding to x_j to the result 
			arr[index(j, column, dim)] = (colLcm == 0 ? 1 : colLcm);
			for (int k = 0; k < std::min(row + 1, n); ++k) {
				if (phi(k, j) != 0) {
					arr[index(indices[k], column, dim)] = phi(k, j) * colLcm / phi(k, indices[k]) * (-1);
				}
			}
			++column;
		} else {
			indices[row++] = j;
		}
	}
	return Matrix<R>(arr, m, dim);
}

/**
 * @brief Complete the span of the column vectors of this matrix to the span of the column vectors of rhs
 * 
 * @tparam R The type of the underlying euclidean ring
 * @param rhs The matrix to complete to
 * 
 * @return A matrix whose columns form a basis for the complement of the span of the columns of this in the span of the columns of rhs
 */
template <typename R>
Matrix<R> Matrix<R>::complete(const Matrix<R>& rhs) const {
	if (*this == 0) return rhs;
    
	int fRank = getRank();
	int sRank = rhs.getRank();
    int rk = fRank;
	int newM = sRank - fRank;
    int cnt = 0;

    R* arr = util::allocate<R>(n * newM);
    for (int j = 0; j < rhs.getM(); ++j) {
        // Add j-th vector of rhs and check if the rank has increased
        Matrix<R> span = mergeHorizontal(rhs.getColumns(0, j));
        int newRk = span.getRank();
        if (newRk > rk) {
            // Rank has increased: Add j-th vector to result
            for (int i = 0; i < n; ++i) {
                arr[index(i, cnt, newM)] = rhs(i, j);
            }
            ++cnt;
			rk = newRk;
        }
		if (cnt == newM) break;
    }
	return Matrix<R>(arr, n, newM);
}

/**
 * @brief Merge this matrix with rhs horizontally
 * 
 * @tparam R The underlying type
 * @param rhs The matrix to be merged with
 * @return The merged matrix
 */
template <typename R>
Matrix<R> Matrix<R>::mergeHorizontal(const Matrix<R>& rhs) const {
	int newM = m + rhs.getM();
	R* arr = util::allocate<R>(n * newM);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			arr[index(i, j, newM)] = (*this)(i, j);
		}
		for (int j = 0; j < rhs.getM(); ++j) {
			arr[index(i, m + j, newM)] = rhs(i, j);
		}
	}
	return Matrix<R>(arr, n, newM);
}

/**
 * @brief Calculate the characteristic polynomial of this matrix
 * 
 * @tparam R The underlying euclidean ring with a static function characteristic()
 * 
 * @note The algorithm used depends on the characteristic of the ring
 * 
 * @return The characteristic polynomial 
 */
template <typename R>
Polynomial<R> Matrix<R>::charPoly() const {
	if (n != m) throw std::invalid_argument("Characteristic polynomial undefined for non-square matrix.");
	if (R::characteristic() == 0) {
		return faddeevLeVerrier();
	} else {
		return charPolyNaive();
	}
}

//----------------------------------------------------------------------|
// Getters
//----------------------------------------------------------------------|

template<typename R>
int Matrix<R>::getRank() const {
	if(reduced) return rank;
	// Reduce a copy of matrix to set the rank
	Matrix<R> A(*this);
	A.rref(false);
	return A.getRank();
}

template<typename R>
R Matrix<R>::getDeterminant() const {
	if(reduced) return determinant;
	// Reduce a copy of matrix to set the determinant
	Matrix<R> A(*this);
	A.rref(false);
	return A.getDeterminant();
}

/**
 * @brief Return the i1-th to i2-th row
 * 
 * @tparam R The underlying type
 * @param i1 The lower row index
 * @param i2 The upper row index (inclusive)
 * 
 * @return The i1-th to i2-th row of this matrix
 */
template<typename R>
Matrix<R> Matrix<R>::getRows(int i1, int i2) const {
	int newN = i2 - i1 + 1;
	return Matrix<R>(util::copy<R, R>(entries + i1, newN * m), newN, m);
}

/**
 * @brief Return the j1-th to j2-th column
 * 
 * @tparam R The underlying type
 * @param j1 The lower column index
 * @param j2 The upper column index (inclusive)
 * @return The j1-th to j2-th column of this matrix
 */
template<typename R>
Matrix<R> Matrix<R>::getColumns(int j1, int j2) const {
	int newM = j2 - j1 + 1;
    R* arr = util::allocate<R>(n * newM);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < newM; ++j) {
            arr[index(i, j, newM)] = (*this)(i, j1 + j);
        }
    }
	return Matrix<R>(arr, n, newM);
}

template <typename R>
void swap(Matrix<R>& lhs, Matrix<R>& rhs) {
	std::swap(lhs.reduced, rhs.reduced);
	std::swap(lhs.determinant, rhs.determinant);
	std::swap(lhs.n, rhs.n);
	std::swap(lhs.m, rhs.m);
	std::swap(lhs.entries, rhs.entries);
}

//-------------------------------------------------------------------------------------------------------------|
// Static functions
//-------------------------------------------------------------------------------------------------------------|

/**
 * @brief Parses a string to a matrix
 * 
 * @tparam R The underlying ring
 * @param[in] str The string to be parsed
 * @param[out] res The resulting matrix
 * @param[in] n The number of rows in the output matrix (if -1, assumes square matrix)
 * 
 * @note The string must be a sequence of integers separated by whitespace
 * 
 * @return true if parsing was successful, false otherwise
 */
template<typename R>
bool Matrix<R>::parseToMatrix(const std::string& str, Matrix<R>& res, int n) {
	int m;
	int num = 0;
	std::istringstream is(str);
    std::vector<int> data;
	while (is >> num) data.emplace_back(num);
    
	int numTotal = int(size(data));
	if (n == -1) {
        // Get dimension by assuming m == n
		n = int(sqrt(numTotal));
		m = n;
	} else {
		m = numTotal / n;
	}
	if (n * m != numTotal) return false;

    R* values = util::allocate<R>(n * m);

    int cnt = 0;
    for (int value : data) values[cnt++] = value;

    res = Matrix<R>(values, n, m);
    return true;
}