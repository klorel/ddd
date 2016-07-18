#pragma once

#include "ComplexPolynomial.h"

class Constraint
{
public:
	Constraint();
	~Constraint();
public:
	ComplexNumber const & lb()const;
	ComplexNumber const & ub()const;
	ComplexPolynomial const & f()const;

	ComplexNumber & lb();
	ComplexNumber & ub();
	ComplexPolynomial & f();

	void print(std::ostream &)const;
	void print(std::ostream &, Problem const & rhs)const;
private:
	ComplexPolynomial _f;
	ComplexNumberPtr _lb;
	ComplexNumberPtr _ub;
};

