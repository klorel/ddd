#pragma once
#include "RealMonomial.h"

class ComplexMonomial {
public:
	friend ComplexMonomialPtr operator+(ComplexMonomial const & lhs, ComplexMonomial const & rhs);
	enum Z_ORDER {
		Z,
		ZH,
	};
public:
	RealMonomial & zH();
	RealMonomial const & zH()const;

	RealMonomial & z();
	RealMonomial const & z()const;

	bool zero() const;

	std::ostream & print(std::ostream & stream)const;
	std::ostream & print(std::ostream & stream, Problem const &)const;
	ComplexMonomialPtr conjugate()const;
	
private:
	std::tuple<RealMonomialPtr, RealMonomialPtr> _non_zero;
public:
	ComplexMonomial();
public:
	static ComplexMonomialPtr ZeroPtr;
	static ComplexMonomialPtr Build(PosInt);
	static ComplexMonomialPtr BuildH(PosInt);
};