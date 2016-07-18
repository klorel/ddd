#include "ComplexPolynomial.h"


std::ostream & operator<<(std::ostream & stream, ComplexPolynomial const & rhs) {
	rhs.print(stream);
	return stream;
}




ComplexPolynomial operator+(ComplexPolynomial const & rhs) {
	ComplexPolynomial result;
	result.insert(rhs, +1);
	return result;
}
ComplexPolynomial operator-(ComplexPolynomial const & rhs) {
	ComplexPolynomial result;
	result.insert(rhs, -1);
	return result;
}

ComplexPolynomial operator+(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs) {
	ComplexPolynomial result;
	result.insert(lhs, +1);
	result.insert(rhs, +1);
	return result;
}
ComplexPolynomial operator-(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs) {
	ComplexPolynomial result;
	result.insert(lhs, +1);
	result.insert(rhs, -1);
	return result;
}
ComplexPolynomial operator*(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs) {
	ComplexPolynomial result;
	for (auto const & lhs_term : lhs.terms()) {
		for (auto const & rhs_term : rhs.terms()) {
			ComplexMonomial const & lhs_monomial(*lhs_term.first);
			ComplexMonomial const & rhs_monomial(*rhs_term.first);
			// z alpha zH alphaH z beta zH betaH = z alpha+beta zH alphaH+betaH
			// can be largely optimize for dense 
			result.terms()[lhs_monomial + rhs_monomial] = lhs_term.second*rhs_term.second;
		}
	}

	return result;
}
ComplexPolynomial operator/(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs) {
	ComplexPolynomial result;
	if (!rhs.isConstant()) {
		throw std::invalid_argument("in operator/(ComplexPolynomial const &lhs, ComplexPolynomial const &rhs) rhs has a non zero degree");
	}
	else {
		ComplexNumber const constant(rhs.constant());
		if (isZero(constant)) {
			throw std::invalid_argument("in operator/(ComplexPolynomial const &lhs, ComplexPolynomial const &rhs): rhs is zero");
		}
		else {
			result.insert(lhs, ComplexNumber(1, 0) / constant);
		}
	}
	return result;
}


ComplexPolynomial ComplexPolynomial::H()const {
	return conjugate();
}
ComplexPolynomial ComplexPolynomial::conjugate()const {
	ComplexPolynomial result;
	for (ComplexTerms::const_reverse_iterator rit(_terms->rbegin()); rit != _terms->rend(); ++rit) {
		(*result._terms)[rit->first->conjugate()] = std::conj(rit->second);
	}
	return result;
}

ComplexPolynomial::ComplexPolynomial() :_terms(new ComplexTerms) {

}
ComplexPolynomial::ComplexPolynomial(Number value) : _terms(new ComplexTerms) {
	terms()[ComplexMonomial::ZeroPtr] = value;
}
ComplexPolynomial::ComplexPolynomial(Number real, Number imag) : _terms(new ComplexTerms) {
	terms()[ComplexMonomial::ZeroPtr] = ComplexNumber(real, imag);
}
ComplexPolynomial::ComplexPolynomial(ComplexNumber const & value) : _terms(new ComplexTerms) {
	terms()[ComplexMonomial::ZeroPtr] = value;
}
ComplexPolynomial ComplexPolynomial::i() {
	ComplexPolynomial result(0, 1);
	return result;

}
ComplexPolynomial ComplexPolynomial::Build(PosInt id) {
	ComplexPolynomial result;
	result.terms()[ComplexMonomial::Build(id)] = 1;
	return result;
}
ComplexPolynomial ComplexPolynomial::Build(PosInt id, Number value) {
	ComplexPolynomial result;
	result.terms()[ComplexMonomial::Build(id)] = value;
	return result;
}
ComplexPolynomial ComplexPolynomial::Build(PosInt id, ComplexNumber const & value) {
	ComplexPolynomial result;
	result.terms()[ComplexMonomial::Build(id)] = value;
	return result;
}

ComplexPolynomial ComplexPolynomial::BuildH(PosInt id) {
	ComplexPolynomial result;
	result.terms()[ComplexMonomial::BuildH(id)] = 1;
	return result;
}
ComplexPolynomial ComplexPolynomial::BuildH(PosInt id, Number value) {
	ComplexPolynomial result;
	result.terms()[ComplexMonomial::BuildH(id)] = value;
	return result;
}
ComplexPolynomial ComplexPolynomial::BuildH(PosInt id, ComplexNumber const & value) {
	ComplexPolynomial result;
	result.terms()[ComplexMonomial::BuildH(id)] = value;
	return result;
}
std::vector<ComplexPolynomial> ComplexPolynomial::BuildVector(PosInt size) {
	std::vector<ComplexPolynomial> result(size);
	for (PosInt i(0); i < size; ++i)
		result[i] = Build(i);
	return result;
}
std::vector<ComplexPolynomial> ComplexPolynomial::BuildVectorH(PosInt size) {
	std::vector<ComplexPolynomial> result(size);
	for (PosInt i(0); i < size; ++i)
		result[i] = BuildH(i);
	return result;
}

ComplexTerms & ComplexPolynomial::terms() {
	return *_terms;
}

ComplexTerms const & ComplexPolynomial::terms() const {
	return *_terms;
}

bool ComplexPolynomial::isConstant()const {
	return _terms->size() == 1 && _terms->begin()->first->zero();
}
ComplexNumber ComplexPolynomial::constant()const {
	return _terms->begin()->second;
}

void ComplexPolynomial::insert(ComplexPolynomial const & rhs, ComplexNumber factor) {
	for (auto const & term : rhs.terms()) {
		terms()[term.first] += factor*term.second;
	}
}
std::ostream & ComplexPolynomial::print(std::ostream & stream) const {
	for (auto const & term : terms()) {
		stream << format(term.second);
		stream << term.first;
	}
	return stream;
}
std::ostream & ComplexPolynomial::print(std::ostream & stream, Problem const & rhs) const {
	for (auto const & term : terms()) {
		stream << format(term.second);
		term.first->print(stream, rhs);
	}
	return stream;
}

void ComplexPolynomial::operator+=(ComplexPolynomial const & rhs) {
	ComplexPolynomial result = (*this) + rhs;
	_terms = result._terms;
}

void ComplexPolynomial::operator-=(ComplexPolynomial const & rhs) {
	ComplexPolynomial result = (*this) - rhs;
	_terms = result._terms;
}

void ComplexPolynomial::operator*=(ComplexPolynomial const & rhs) {
	ComplexPolynomial result = (*this) * rhs;
	_terms = result._terms;
}

void ComplexPolynomial::operator/=(ComplexPolynomial const & rhs) {
	ComplexPolynomial result = (*this) / rhs;
	_terms = result._terms;
}