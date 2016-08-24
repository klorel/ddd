#pragma once

#include "PolynomialOptimizationProblem.h"

enum AvailablePopInstances {
	//PROBLEM_2_2,
	PROBLEM_2_6,
	PROBLEM_2_9,
	PROBLEM_3_3,
	PROBLEM_3_4,
	PROBLEM_4_6,
	PROBLEM_4_7,
};

template<AvailablePopInstances __pop_instance__> inline void GetInstance(PolynomialOptimizationProblem & pop) {

}


template<> inline void GetInstance<PROBLEM_2_6>(PolynomialOptimizationProblem & pop) {

	NumberVector c = { 48, 42, 48, 45, 44, 41, 47, 42, 45, 4 };
	NumberVector b = { -4,22,-6,-23,-12 };

	NumberVector A = {
		-2,-6,-1, 0,-3,-3,-2,-6,-2,-2,
		 6,-5, 8,-3, 0, 1, 3, 8, 9,-3,
		-5, 6, 5, 3, 8,-8, 9, 2, 0,-9,
		 9, 5, 0,-9, 1,-8, 3,-9,-9,-3,
		-8, 7,-4,-5,-9, 1,-7,-1, 3,-2 };

	Number Q = 100;
	pop.clear();
	int const n(10);
	IndexedPool const & x = pop.newvarpool("x", n);
	for (int i(0); i < n; ++i) {
		std::cout << "nb Term : "<<pop.minimize().terms().size() << std::endl;
		pop.minimize() += c[i] * x(i) -0.5*Q*x(i)*x(i);
		pop.add(x(i)*x(i) <= 1);
		//pop.add(x(i)*x(i) >= 0);
	}
	for (int i(0); i < (int)b.size(); ++i) {
		ComplexPolynomial f;
		for (int j(0); j < n; ++j) {
			f += A[i*n + j] * x(j);
		}
		//pop.add(f <= b[i]);
	}
}