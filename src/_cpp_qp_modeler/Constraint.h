#pragma once

#include "FunctionReal.h"

class Constraint
{
public:
	Constraint();
	~Constraint();
public:
	Number const & lb()const;
	Number const & ub()const;
	FunctionReal const & f()const;

	Number & lb();
	Number & ub();
	FunctionReal & f();

	void print(std::ostream &)const;
	void print(std::ostream &, PolynomialOptimizationProblem const &)const;

	Sense sense()const;
private:
	FunctionReal _f;
	NumberPtr _lb;
	NumberPtr _ub;
};

