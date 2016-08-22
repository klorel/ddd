#pragma once
#include "common.h"

class RealMonomial {
public:
	friend RealMonomialPtr operator+(RealMonomial const & lhs, RealMonomial const & rhs);
public:
	PosInt2PosInt & alpha();
	PosInt2PosInt const & alpha()const;

	RealMonomial();
	PosInt degree()const;
	PosInt &degree();

	std::ostream & print(std::ostream & stream, bool const & isH)const;
	std::ostream & print(std::ostream & stream, PolynomialOptimizationProblem const & rhs, bool const & isH)const;
	
private:
	PosInt2PosInt _alpha;
	PosInt _degree;
private:
};


//class RealMonomial;
//class ComplexMonomial;
//typedef std::shared_ptr<ComplexMonomial> ComplexMonomialPtr;
//typedef std::shared_ptr<RealMonomial> RealMonomialPtr;
//std::ostream & operator<<(std::ostream & stream, ComplexMonomial const & rhs);
//
//#define __DEGUB_ORDERING__ 0
//RealMonomialPtr operator+(RealMonomial const & lhs, RealMonomial const & rhs);
//
//class RealMonomial {
//public:
//	friend RealMonomialPtr operator+(RealMonomial const & lhs, RealMonomial const & rhs);
//public:
//	PosInt2PosInt & alpha() { return _alpha; }
//	PosInt2PosInt const & alpha()const { return _alpha; }
//	
//	RealMonomial() :_alpha(), _degree(0) {
//
//	}
//	PosInt degree()const { return _degree; }
//	PosInt &degree() { return _degree; }
//
//	std::ostream & print(std::ostream & stream, bool const & isH)const {
//		for (auto const & alpha : alpha()) {
//			stream << "z" << (isH ? "H" : "") << "[" << alpha.first << "]";
//			if (alpha.second > 1)
//				stream << "^" << alpha.second;
//		}
//		return stream;
//	}
//private:
//	PosInt2PosInt _alpha;
//	PosInt _degree;
//private:
//};
//inline std::ostream & operator<<(std::ostream & stream, RealMonomial const & rhs) {
//	rhs.print(stream, true);
//	return stream;
//}
//inline bool operator<(RealMonomial const &  lhs, RealMonomial const & rhs) {
//	bool result;
//	if (lhs.alpha().empty())
//		result = false;
//	else if (rhs.alpha().empty())
//		result = true;
//	else if (lhs.degree() == rhs.degree())
//		result = lhs.alpha() < rhs.alpha();
//	else
//		result = lhs.degree() < rhs.degree();
//#if __DEGUB_ORDERING__
//	std::cout << "RealMonomial "<<lhs << (result ? "<" : "=>") << rhs << std::endl;
//#endif
//	return result;
//}
//
//inline bool operator>(RealMonomial const &  lhs, RealMonomial const & rhs) {
//	bool result;
//	if (rhs.alpha().empty())
//		result = false;
//	else if (lhs.alpha().empty())
//		result = true;
//	else if (lhs.degree() == rhs.degree())
//		result = lhs.alpha() < rhs.alpha();
//	else
//		result = lhs.degree() < rhs.degree();
//#if __DEGUB_ORDERING__
//	std::cout << "RealMonomial " << lhs << (result ? ">" : "<=") << rhs << std::endl;
//#endif
//	return result;
//}
//
//inline bool operator==(RealMonomial const &  lhs, RealMonomial const & rhs) {
//	bool result;
//	if (lhs.degree() == rhs.degree())
//		result = lhs.alpha() == rhs.alpha();
//	else
//		result = false;
//#if __DEGUB_ORDERING__
//	std::cout << "RealMonomial " << lhs << (result ? "==" : "!=") << rhs << std::endl;
//#endif
//	return result;
//}
