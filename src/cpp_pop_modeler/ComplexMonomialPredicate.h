#pragma once
#include "ComplexMonomial.h"

class ComplexMonomialPredicate {
public:
	bool operator()(ComplexMonomial const & lhs, ComplexMonomial const & rhs)const;
	bool operator()(ComplexMonomialPtr const & lhs, ComplexMonomialPtr const & rhs) const {
		return operator()(*lhs, *rhs);
	}
};