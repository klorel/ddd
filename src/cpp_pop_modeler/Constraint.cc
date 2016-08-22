#include "Constraint.h"

#include "PolynomialOptimizationProblem.h"

Constraint::Constraint() :_lb(new ComplexNumber(negInfinity())), _ub(new ComplexNumber(posInfinity())), _f()
{
}


Constraint::~Constraint()
{
}

ComplexNumber const & Constraint::lb()const {
	return *_lb;
}
ComplexNumber const & Constraint::ub()const {
	return *_ub;
}
ComplexPolynomial const & Constraint::f()const {
	return _f;
}

ComplexNumber & Constraint::lb() {
	return *_lb;
}
ComplexNumber & Constraint::ub() {
	return *_ub;
}
ComplexPolynomial & Constraint::f() {
	return _f;
}
void Constraint::print(std::ostream & stream)const {
	stream << lb() << " <= " << f() << " <= " << ub();
}
void Constraint::print(std::ostream & stream, PolynomialOptimizationProblem const & rhs)const {
	if (lb() == ub()) {
		stream << lb() << "  = ";
		f().print(stream, rhs);
	}
	else {
		stream << lb() << " <= ";
		f().print(stream, rhs);
		stream << " <= " << ub();
	}
}
std::ostream & operator<<(std::ostream & lhs, Constraint const & rhs) {
	rhs.print(lhs);
	return lhs;
} 

Constraint operator<=(ComplexPolynomial const & lhs, ComplexNumber const & rhs) {
	Constraint result;
	result.f() = lhs;
	result.ub() = rhs;
	return result;

}
Constraint operator==(ComplexPolynomial const & lhs, ComplexNumber const & rhs) {
	Constraint result;
	result.f() = lhs;
	result.ub() = rhs;
	result.lb() = rhs;
	return result;

}
Constraint operator>=(ComplexPolynomial const & lhs, ComplexNumber const & rhs) {
	Constraint result;
	result.f() = lhs;
	result.lb() = rhs;
	return result;
}

Constraint operator<=(ComplexNumber const & lhs, ComplexPolynomial const & rhs) {
	Constraint result;
	result.f() = rhs;
	result.lb() = lhs;
	return result;
}

Constraint operator==(ComplexNumber const & lhs, ComplexPolynomial const & rhs) {
	Constraint result;
	result.f() = rhs;
	result.ub() = lhs;
	result.lb() = lhs;
	return result;
}

Constraint operator>=(ComplexNumber const & lhs, ComplexPolynomial const & rhs) {
	Constraint result;
	result.f() = rhs;
	result.ub() = lhs;
	return result;
}