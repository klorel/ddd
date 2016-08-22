#include "FunctionComplex.h"

FunctionComplex::FunctionComplex() {
}
FunctionComplex::FunctionComplex(Number v1, Number v2) :_real(v1), _imag(v2) {

}
FunctionComplex::~FunctionComplex() {
}


FunctionReal const & FunctionComplex::real()const {
	return _real;
}
FunctionReal & FunctionComplex::real() {
	return _real;
}

FunctionReal const & FunctionComplex::imag()const {
	return _imag;
}
FunctionReal & FunctionComplex::imag() {
	return _imag;
}

FunctionComplex FunctionComplex::clone()const {
	FunctionComplex result;
	result.real() = real().clone();
	result.imag() = imag().clone();
	return result;
}
FunctionComplex FunctionComplex::conjugate()const {
	FunctionComplex result;
	result.real() = real().clone();
	result.imag() = -imag().clone();
	return result;
}
bool FunctionComplex::isReal()const {
	return _imag.isNull();
}
bool FunctionComplex::isImag()const {
	return _real.isNull();
}
bool FunctionComplex::isNull()const {
	return isReal() && isImag();
}

bool FunctionComplex::isConstant()const {
	return real().isConstant() && imag().isConstant();
}

bool FunctionComplex::isLinear()const {
	return real().isLinear() && imag().isLinear();
}

void FunctionComplex::print(std::ostream & stream, PolynomialOptimizationProblem const & problem)const {
	if (!isNull()) {
		if (!real().isNull()) {
			real().print(stream, problem);
		}
		if (!imag().isNull()) {
			imag().print(stream << "+i*(", problem);
			stream << ")";
		}
	}
	else {
		stream << "0";
	}
}
void FunctionComplex::print(std::ostream & stream)const {
	if (!isNull()) {
		if (!real().isNull()) {
			real().print(stream);
		}
		if (!imag().isNull()) {
			imag().print(stream << "+i*(");
			stream << ")";
		}
	}
	else {
		stream << "0";
	}
}


void FunctionComplex::operator+=(FunctionComplex const & rhs) {
	real() += rhs.real().clone();
	imag() += rhs.imag().clone();
}
void FunctionComplex::operator-=(FunctionComplex const & rhs) {
	real() -= rhs.real().clone();
	imag() -= rhs.imag().clone();

}
void FunctionComplex::operator*=(FunctionComplex const & rhs) {
	FunctionReal r = real().clone()*rhs.real().clone() - imag().clone()*rhs.imag().clone();
	FunctionReal i = imag().clone() = real()*rhs.imag().clone() + imag().clone()*rhs.real().clone();
	real() = r;
	imag() = i;

}
FunctionComplex FunctionComplex::module2()const {
	return clone()*conjugate();
}
void FunctionComplex::operator/=(FunctionComplex const & rhs) {
	if (rhs.isNull() || !rhs.isConstant()) {
		throw std::invalid_argument("in FunctionComplex::operator/=");
	}
	else {
		FunctionComplex m2 = rhs.module2();
		if (!m2.isReal() && !m2.isConstant())
			throw std::invalid_argument("in FunctionComplex::operator/=");
		Number m2value = m2.real().constant();
		real() /= m2value;
		imag() /= m2value;

		operator*=( rhs.conjugate() );
	}
}

