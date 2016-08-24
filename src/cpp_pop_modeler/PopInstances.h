#pragma once

#include "PolynomialOptimizationProblem.h"

enum AvailablePopInstances {
	//PROBLEM_2_2,
	EXAMPLE_2,
	PROBLEM_2_6,
	PROBLEM_2_9,
	PROBLEM_3_3,
	PROBLEM_3_4,
	PROBLEM_4_6,
	PROBLEM_4_7,
};

template<AvailablePopInstances __pop_instance__> inline void GetInstance(PolynomialOptimizationProblem & pop) {

}

template<> inline void GetInstance<EXAMPLE_2>(PolynomialOptimizationProblem & pop) {
	IndexedPool const & x = pop.newvarpool("x", 2);
	ComplexPolynomial term1 = (1 + x(0)*x(0));
	ComplexPolynomial term2 = (1 + x(1)*x(1));
	ComplexPolynomial term3 = (1 + x(0) + x(1));
	pop.minimize() = -2 * term3*term3;
	pop.minimize() = term1*term1 + term2*term2 - 2 * term3*term3;

	double const x_opt(1.3247);
	double const x_lb(1e-2 + 1.3247);
	pop.add(x(0)*x(0) >= x_lb*x_lb);
	pop.add(x(1)*x(1) >= x_lb*x_lb);
}

template<> inline void GetInstance<PROBLEM_2_6>(PolynomialOptimizationProblem & pop) {
	pop.clear();

	NumberVector c = { 48, 42, 48, 45, 44, 41, 47, 42, 45, 4 };
	NumberVector b = { -4,22,-6,-23,-12 };

	NumberVector A = {
		-2,-6,-1, 0,-3,-3,-2,-6,-2,-2,
		 6,-5, 8,-3, 0, 1, 3, 8, 9,-3,
		-5, 6, 5, 3, 8,-8, 9, 2, 0,-9,
		 9, 5, 0,-9, 1,-8, 3,-9,-9,-3,
		-8, 7,-4,-5,-9, 1,-7,-1, 3,-2 };

	Number Q = 100;
	int const n(10);
	IndexedPool const & x = pop.newvarpool("x", n);
	for (int i(0); i < n; ++i) {
		pop.minimize() += c[i] * x(i) - 0.5*Q*x(i)*x(i);
		pop.add(x(i)*x(i) <= 1);
		pop.add(x(i)*x(i) >= 0);
	}
	for (int i(0); i < (int)b.size(); ++i) {
		ComplexPolynomial f;
		for (int j(0); j < n; ++j) {
			f += A[i*n + j] * x(j);
		}
		pop.add(f <= b[i]);
	}

}
template<> inline void GetInstance<PROBLEM_4_6>(PolynomialOptimizationProblem & pop) {
	pop.clear();
	int const n(2);
	IndexedPool const & x = pop.newvarpool("x", n);
	pop.minimize() += -x(0) - x(1);

	pop.add(0 <= -x(1) + 2 * ComplexPolynomial::Pow(0, 4) - 8 * ComplexPolynomial::Pow(0, 3) + 8 * ComplexPolynomial::Pow(0, 2) + 2);
	pop.add(0 <= -x(1) + 4 * ComplexPolynomial::Pow(0, 4) - 32 * ComplexPolynomial::Pow(0, 3) + 88 * ComplexPolynomial::Pow(0, 2) - 96 * x(0) + 36);

	pop.add(x(0) <= 3);
	pop.add(x(0) >= 0);

	pop.add(x(1) <= 4);
	pop.add(x(1) >= 0);


}
template<> inline void GetInstance<PROBLEM_4_7>(PolynomialOptimizationProblem & pop) {
	pop.clear();
	int const n(2);
	IndexedPool const & x = pop.newvarpool("x", n);
	pop.minimize() += -12*x(0) - 7*x(1)+ ComplexPolynomial::Pow(1, 2);

	pop.add(0 <= -2* ComplexPolynomial::Pow(0, 4)+2-x(1));
	pop.add(0 >= -2 * ComplexPolynomial::Pow(0, 4) + 2 - x(1));

	pop.add(x(0) <= 2);
	pop.add(x(0) >= 0);

	pop.add(x(1) <= 3);
	pop.add(x(1) >= 0);


}
template<> inline void GetInstance<PROBLEM_3_3>(PolynomialOptimizationProblem & pop) {
	pop.clear();
	int const n(6);
	IndexedPool const & x = pop.newvarpool("x", n);
	pop.add(x(0) >= 0); pop.add(x(0)*x(0) <= 1000);
	pop.add(x(1) >= 0); pop.add(x(1)*x(1) <= 1000);
	pop.add(x(2) >= 1); pop.add(x(2)*x(2) <= 25);
	pop.add(x(3) >= 0); pop.add(x(3)*x(3) <= 36);
	pop.add(x(4) >= 1); pop.add(x(4)*x(4) <= 25);
	pop.add(x(5) >= 0); pop.add(x(5)*x(5) <= 100);

	ComplexPolynomial term_1 = x(0) - 2;
	ComplexPolynomial term_2 = x(1) - 2;
	ComplexPolynomial term_3 = x(2) - 1;
	ComplexPolynomial term_4 = x(3) - 4;
	ComplexPolynomial term_5 = x(4) - 1;
	ComplexPolynomial term_6 = x(5) - 4;

	pop.minimize() += -25 * term_1*term_1;
	pop.minimize() += - term_2*term_2;
	pop.minimize() += - term_3*term_3;
	pop.minimize() += - term_4*term_4;
	pop.minimize() += - term_5*term_5;
	pop.minimize() += - term_6*term_6;

	pop.add(4 <= (x(2) - 3)*(x(2) - 3) + x(3));
	pop.add(4 <= (x(4) - 3)*(x(4) - 3) + x(5));

	pop.add(x(0) - 3 * x(1) <= 2);
	pop.add(-x(0) + x(1) <= 2);
	pop.add(x(0) + x(1) <= 6);
	pop.add(x(0) + x(1) >= 2);
}