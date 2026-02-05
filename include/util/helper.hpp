#pragma once
#include <string>

namespace util {

	static int numAlloc = 0;
	static int numDealloc = 0;
	static int sizeAlloc = 0;

	/**
	* @brief Greatest common divisor of two numbers
	*
	* @tparam R A class representing any euclidean ring
	*
	* @param a, b Two numbers to compute the gcd of
	*
	* @return the gcd of a and b
	*/
	template <typename R>
	inline R gcd(const R& a, const R& b) {
		R result = (b == 0) ? a : gcd(b, a % b);
		return result;
	}
	
	/**
	* @brief A custom wrapper function of allocation using new
	*
	* @tparam T The type of the array to allocate
	*
	* @param k The size of the array to allocate
	*
	* @return a pointer to the allocated array
	*/
	template <typename T>
	T* allocate(int k) {
		if (k > 0) {
			++numAlloc;
			sizeAlloc += sizeof(T) * k;
			return new T[k];
		}
		return nullptr;
	}
	
	/**
	 * @brief A custom wrapper function of copying an array of one type to another
	 * 
	 * @tparam T_IN The input type of the array to copy
	 * @tparam T_OUT The output type of the array to copy
	 * 
	 * @param arr The input array to copy
	 * @param k The size of the input array to copy
	 * 
	 * @return a pointer to the copied array of type T_OUT
	 */
	template <typename T_IN, typename T_OUT>
	T_OUT* copy(T_IN* arr, int k) {
		T_OUT* res = allocate<T_OUT>(k);
		std::copy(arr, arr + k, res);
		return res;
	}

	/**
	 * @brief A function allocating an array of zeroes
	 * 
	 * @tparam R A class representing any ring
	 * 
	 * @param k The size of the array to allocate
	 * 
	 * @return a pointer to the allocated array of zeroes 
	 */
	template <typename R>
	R* zeroes(int k) {
		R* res = allocate<R>(k);
		for (int j = 0; j < k; ++j) res[j] = 0;
		return res;
	}

	/**
	* @brief A wrapper function of delete[]
	*
	* @tparam T The type of the array to deallocate
	*
	* @param ptr A pointer to the array to deallocate
	*/
	template <typename T>
	void deallocate(T* ptr) {
		if (ptr != nullptr) {
			++numDealloc;
			delete[] ptr;
		}
	}

	/**
	 * @brief A function to assert equality of two values and print an error message if they are not
	 * 
	 * @tparam T_VAL The type of the value to be tested
	 * @tparam T_GTRUTH The type of the expected value
	 * @tparam T_OBJ The type of the object being tested
	 * 
	 * @param val The value to be tested
	 * @param gTruth The expected value
	 * @param obj The object being tested
	 * @param msg A message to be printed in case of failure
	 * 
	 * @return true if the values are equal, false otherwise
	 */
    template<typename T_VAL, typename T_GTRUTH, typename T_OBJ>
    bool assertEq(T_VAL val, T_GTRUTH gTruth, T_OBJ obj, std::string msg) {
		if (val != gTruth) {
			std::cerr << "Error testing " << msg << ":\n" << val << "\nnot equal to\n" << gTruth <<"\nOccured using\n" << obj << "\n";
			return EXIT_FAILURE;
		} 
		return EXIT_SUCCESS;
	}
}