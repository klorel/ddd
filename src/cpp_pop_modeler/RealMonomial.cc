#include "RealMonomial.h"

#include "PolynomialOptimizationProblem.h"

PosInt2PosInt & RealMonomial::alpha() {
	return _alpha;
}
PosInt2PosInt const & RealMonomial::alpha()const {
	return _alpha;
}

RealMonomial::RealMonomial() :_alpha(), _degree(0) {

}
PosInt RealMonomial::degree()const {
	return _degree;
}
PosInt & RealMonomial::degree() {
	return _degree;
}

std::ostream & RealMonomial::print(std::ostream & stream, PolynomialOptimizationProblem const & rhs,  bool const & isH)const {
	for (auto const & alpha : alpha()) {
		stream << rhs.name(alpha.first) << (isH ? "H" : "");
		if (alpha.second > 1)
			stream << "^" << alpha.second;
	}
	return stream;
}
std::ostream & RealMonomial::print(std::ostream & stream, bool const & isH)const {
	for (auto const & alpha : alpha()) {
		stream << "z" << "[" << alpha.first << "]" << (isH ? "H" : "");
		if (alpha.second > 1)
			stream << "^" << alpha.second;
	}
	return stream;
}


std::ostream & operator<<(std::ostream & stream, RealMonomial const & rhs) {
	rhs.print(stream, true);
	return stream;
}
bool operator<(RealMonomial const &  lhs, RealMonomial const & rhs) {
	bool result;
	if (lhs.alpha().empty())
		result = false;
	else if (rhs.alpha().empty())
		result = true;
	else if (lhs.degree() == rhs.degree())
		result = lhs.alpha() < rhs.alpha();
	else
		result = lhs.degree() < rhs.degree();
#if __DEGUB_ORDERING__
	std::cout << "RealMonomial " << lhs << (result ? "<" : "=>") << rhs << std::endl;
#endif
	return result;
}

bool operator>(RealMonomial const &  lhs, RealMonomial const & rhs) {
	bool result;
	if (rhs.alpha().empty())
		result = false;
	else if (lhs.alpha().empty())
		result = true;
	else if (lhs.degree() == rhs.degree())
		result = lhs.alpha() < rhs.alpha();
	else
		result = lhs.degree() < rhs.degree();
#if __DEGUB_ORDERING__
	std::cout << "RealMonomial " << lhs << (result ? ">" : "<=") << rhs << std::endl;
#endif
	return result;
}

bool operator==(RealMonomial const &  lhs, RealMonomial const & rhs) {
	bool result;
	if (lhs.degree() == rhs.degree())
		result = lhs.alpha() == rhs.alpha();
	else
		result = false;
#if __DEGUB_ORDERING__
	std::cout << "RealMonomial " << lhs << (result ? "==" : "!=") << rhs << std::endl;
#endif
	return result;
}

RealMonomialPtr operator+(RealMonomial const & lhs, RealMonomial const & rhs) {
	// (alpha) + (beta) = (alpha+beta)
	RealMonomialPtr ptr(new RealMonomial);
	RealMonomial & result(*ptr);
	for (auto const & term : lhs.alpha())
		result.alpha()[term.first] += term.second;
	for (auto const & term : rhs.alpha())
		result.alpha()[term.first] += term.second;
	result.degree() = lhs.degree() + rhs.degree();
	return ptr;
}