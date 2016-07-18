#include "ComplexMonomial.h"
#include "Problem.h"

ComplexMonomialPtr  ComplexMonomial::ZeroPtr = ComplexMonomialPtr(new ComplexMonomial);


ComplexMonomialPtr ComplexMonomial::Build(PosInt id) {
	ComplexMonomialPtr result(new ComplexMonomial);
	result->z().alpha()[id] = 1;
	result->z().degree() = 1;
	return result;
}
ComplexMonomialPtr ComplexMonomial::BuildH(PosInt id) {
	ComplexMonomialPtr result(new ComplexMonomial);
	result->zH().alpha()[id] = 1;
	result->zH().degree() = 1;
	return result;
}

RealMonomial & ComplexMonomial::zH() {
	return *std::get<ZH>(_non_zero);
}
RealMonomial const & ComplexMonomial::zH()const {
	return *std::get<ZH>(_non_zero);
}

RealMonomial & ComplexMonomial::z() {
	return *std::get<Z>(_non_zero);
}
RealMonomial const & ComplexMonomial::z()const {
	return *std::get<Z>(_non_zero);
}

bool ComplexMonomial::zero() const {
	return zH().degree() == 0 && z().degree() == 0;
}

std::ostream & ComplexMonomial::print(std::ostream & stream)const {
	zH().print(stream, true);
	z().print(stream, false);

	return stream;
}
std::ostream & ComplexMonomial::print(std::ostream & stream, Problem const & rhs)const {
	zH().print(stream, rhs, true);
	z().print(stream, rhs, false);

	return stream;
}
ComplexMonomialPtr ComplexMonomial::conjugate()const {
	ComplexMonomialPtr result(new ComplexMonomial);
	result->_non_zero = _non_zero;
	std::swap(std::get<Z>(result->_non_zero), std::get<ZH>(result->_non_zero));
	return result;
}

ComplexMonomial::ComplexMonomial() {
	std::get<ZH>(_non_zero) = RealMonomialPtr(new RealMonomial);
	std::get<Z>(_non_zero) = RealMonomialPtr(new RealMonomial);
}


std::ostream & operator<<(std::ostream & stream, ComplexMonomial const & rhs) {
	rhs.print(stream);
	return stream;
}

std::ostream & operator<<(std::ostream & stream, ComplexMonomialPtr const & rhs) {
	return stream << *rhs;
}

ComplexMonomialPtr operator+(ComplexMonomial const & lhs, ComplexMonomial const & rhs) {
	// (alpha, alphaH) + (beta, betaH) = (alpha+beta, alphaH+betaH)
	ComplexMonomialPtr ptr(new ComplexMonomial);
	ComplexMonomial & result(*ptr);
	std::get<ComplexMonomial::Z>(result._non_zero) = lhs.z() + rhs.z();
	std::get<ComplexMonomial::ZH>(result._non_zero) = lhs.zH() + rhs.zH();
	return ptr;
}