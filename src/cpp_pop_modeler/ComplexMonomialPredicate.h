#pragma once
#include "ComplexMonomial.h"

class ComplexMonomialPredicate {
public:
	bool operator()(ComplexMonomial const & lhs, ComplexMonomial const & rhs)const {
		//return lhs._degreeToid < rhs._degreeToid;
		return lhs._id2degree < rhs._id2degree;
//		bool result;
//		if (lhs.zero())
//			result = !rhs.zero();
//		else if (rhs.zero())
//			result = false;
//		else if (lhs.zH() == rhs.zH()) {
//			result = lhs.z()< rhs.z();			
//		}
//		else {
//			result = lhs.zH() < rhs.zH();
//		}
//#if __DEGUB_ORDERING__
//		std::cout << "ComplexMonomial " << lhs ;
//		std::cout << (result ? " < " : " => ");
//		std::cout << rhs << std::endl;
//#endif
		//return result;		
	}
	bool operator()(ComplexMonomialPtr const & lhs, ComplexMonomialPtr const & rhs) const {
		return operator()(*lhs, *rhs);
	}
};