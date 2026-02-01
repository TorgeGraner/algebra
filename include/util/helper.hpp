#pragma once
#include <cassert>
#include <string>

namespace util {

	static int numAlloc = 0;
	static int numDealloc = 0;
	static int sizeAlloc = 0;

	/*
	* @brief Greatest common divisor of two things
	* @return the gcd of a and b
	* Probably not the fastest way to do this
	*/
	template <typename R>
	inline R gcd(const R& a, const R& b) {
		R result = (b == 0) ? a : gcd(b, a % b);
		return result;
	}
	
	/*
	* @brief A custom wrapper function of allocation using new
	* @param final If set to true print the number of allocated arrays
	* Used to keep track of the number of allocations for testing if the destructors and constructors of the classes work as intended.
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
	
	template <typename T_IN, typename T_OUT>
	T_OUT* copy(T_IN* arr, int k) {
		T_OUT* res = allocate<T_OUT>(k);
		std::copy(arr, arr + k, res);
		return res;
	}

	template <typename R>
	R* zeroes(int k) {
		R* res = allocate<R>(k);
		for (int j = 0; j < k; ++j) res[j] = 0;
		return res;
	}

	/*
	* @brief A wrapper function of delete[]
	* @param final If set to true, print the number of deallocated arrays
	* Used to keep track of the number of deallocations for testing if the destructors and constructors of the classes work as intended.
	*/
	template <typename T>
	void deallocate(T* ptr) {
		if (ptr != nullptr) {
			++numDealloc;
			delete[] ptr;
		}
	}

	template<typename Base, typename T>
	inline bool instanceof(const T*) {
		return std::is_base_of<Base, T>::value;
	}

    template<typename T_VAL, typename T_GTRUTH, typename T_OBJ>
    bool assertEq(T_VAL val, T_GTRUTH gTruth, T_OBJ obj, std::string msg) {
		if (val != gTruth) {
			std::cerr << "Error testing " << msg << ":\n" << val << "\nnot equal to\n" << gTruth <<"\nOccured using\n" << obj << "\n";
			return EXIT_FAILURE;
		} 
		return EXIT_SUCCESS;
	}
}