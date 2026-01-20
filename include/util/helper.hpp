#pragma once
#include <iostream>
#include <string>

namespace util {
	/*
	* @brief Greatest common divisor of two things
	* @return the gcd of a and b
	* Probably not the fastest way to do this
	*/
	template <typename R>
	inline R gcd(const R& a, const R& b) {
		return (b == 0) ? a : gcd(b, a % b);
	}
	
	/*
	* @brief A custom wrapper function of allocation using new
	* @paramn final If set to true print the number of allocated arrays
	* Used to keep track of the number of allocations for testing if the destructors and constructors of the classes work as intended.
	*/
	template <typename R>
	R* allocate(int k, bool final = false) {
		static int cnt = 0;
		static int cntSz = 0;
		if (!final) {
			++cnt;
			cntSz += sizeof(R) * k;
			return new R[k];
		} else {
			std::cout << "Allocated " << cnt << " arrays with in total " << cntSz << " bytes of " << typeid(R).name() << std::endl;
			return nullptr;
		}
	}

	/*
	* @brief A custom wrapper function of deallocation using delete[]
	* @param final If set to true, print the number of deallocated arrays
	* Used to keep track of the number of deallocations for testing if the destructors and constructors of the classes work as intended.
	*/
	template <typename R>
	void deallocate(R* ptr, bool final = false) {
		static int cnt = 0;
		if (!final) {
			if (ptr != nullptr) {
				++cnt;
				delete[] ptr;
			}
		} else {
			std::cout << "Deallocated " << cnt << " arrays of " << typeid(R).name() << std::endl;
		}
	}
}