#pragma once
#include <iostream>
#include <string>
#include <cassert>

namespace util {

	static int numAlloc = 0;
	static int sizeAlloc = 0;

	static int numDealloc = 0;

	/*
	* @brief Greatest common divisor of two things
	* @return the gcd of a and b
	* Probably not the fastest way to do this
	*/
	template <typename R>
	inline R gcd(const R& a, const R& b) {
		R result = (b == 0) ? a : gcd(b, a % b);
		if (a != 0 || b != 0) assert(result != 0);
		return result;
	}
	
	/*
	* @brief A custom wrapper function of allocation using new
	* @paramn final If set to true print the number of allocated arrays
	* Used to keep track of the number of allocations for testing if the destructors and constructors of the classes work as intended.
	*/
	template <typename R>
	R* allocate(int k) {
		if (k > 0) {
			++numAlloc;
			sizeAlloc += sizeof(R) * k;
			return new R[k];
		}
		return nullptr;
	}

	/*
	* @brief A custom wrapper function of deallocation using delete[]
	* @param final If set to true, print the number of deallocated arrays
	* Used to keep track of the number of deallocations for testing if the destructors and constructors of the classes work as intended.
	*/
	template <typename R>
	void deallocate(R* ptr) {
		if (ptr != nullptr) {
			++numDealloc;
			delete[] ptr;
		}
	}
}