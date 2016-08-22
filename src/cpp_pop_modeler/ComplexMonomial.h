#pragma once
#include "RealMonomial.h"

class ComplexMonomial {
public:
	friend ComplexMonomialPtr operator+(ComplexMonomial const & lhs, ComplexMonomial const & rhs);

public:

	bool zero() const;

	std::ostream & print(std::ostream & stream)const;
	std::ostream & print(std::ostream & stream, PolynomialOptimizationProblem const &)const;

	ComplexMonomialPtr conjugate()const;
	
	int degree()const;
//private:
//	std::tuple<RealMonomialPtr, RealMonomialPtr> _non_zero;
public:
	ComplexMonomial();
public:
	static ComplexMonomialPtr ZeroPtr;
	static ComplexMonomialPtr Build(int);
	static ComplexMonomialPtr BuildH(int);
public:
	// id of a variable -- > degree d, z ^ d if d>0 or z H ^ d if d<0
	Int2Int _id2degree;
	Int2IntSet _degreeToid;
};