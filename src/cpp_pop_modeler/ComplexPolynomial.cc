#include "ComplexPolynomial.h"


std::ostream & operator<<(std::ostream & stream, ComplexPolynomial const & rhs) {
	rhs.print(stream);
	return stream;
}

void ComplexPolynomial::get_all_monomial(ComplexMonomialPtr2Int & output)const {
	for (auto const & term : terms()) {
		output.insert({ term.first, 1 });
	}
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
	//std::cout << "lhs = " << lhs << std::endl;
	//std::cout << "rhs = " << rhs << std::endl;
	for (auto const & lhs_term : lhs.terms()) {
		for (auto const & rhs_term : rhs.terms()) {
			ComplexMonomialPtr key(lhs_term.first + rhs_term.first);
			result.insert(key, lhs_term.second*rhs_term.second, 1);
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
			ComplexNumber inv = ComplexNumber(1.0, 0) / constant;
			result.insert(lhs, inv);
		}
	}
	//std::cout << lhs << " / " << rhs << " = "<< result << std::endl;
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

ComplexPolynomial::ComplexPolynomial() :_terms(new ComplexTerms), _degree(new Degree2Terms){

}
ComplexPolynomial::ComplexPolynomial(Number value) : _terms(new ComplexTerms), _degree(new Degree2Terms) {
	insert(ComplexMonomial::ZeroPtr, value, 1);
}
ComplexPolynomial::ComplexPolynomial(Number real, Number imag) : _terms(new ComplexTerms), _degree(new Degree2Terms) {
	insert(ComplexMonomial::ZeroPtr, ComplexNumber(real, imag), 1);
}
ComplexPolynomial::ComplexPolynomial(ComplexNumber const & value) : _terms(new ComplexTerms), _degree(new Degree2Terms) {
	insert(ComplexMonomial::ZeroPtr, value, 1);
}

ComplexPolynomial::ComplexPolynomial(ComplexMonomialPtr ptr) : _terms(new ComplexTerms), _degree(new Degree2Terms) {
	insert(ptr, 1, 1);
}

ComplexPolynomial ComplexPolynomial::i() {
	ComplexPolynomial result(0, 1);
	return result;

}
ComplexPolynomial ComplexPolynomial::Build(PosInt id) {
	ComplexPolynomial result;
	result.insert(ComplexMonomial::Build(id), 1, 1);
	return result;
}
ComplexPolynomial ComplexPolynomial::Build(PosInt id, Number value) {
	ComplexPolynomial result;
	result.insert(ComplexMonomial::Build(id), value, 1);
	return result;
}
ComplexPolynomial ComplexPolynomial::Build(PosInt id, ComplexNumber const & value) {
	ComplexPolynomial result;
	result.insert(ComplexMonomial::Build(id), value, 1);
	return result;
}

ComplexPolynomial ComplexPolynomial::BuildH(PosInt id) {
	ComplexPolynomial result;
	result.insert(ComplexMonomial::BuildH(id), 1, 1);
	return result;
}
ComplexPolynomial ComplexPolynomial::BuildH(PosInt id, Number value) {
	ComplexPolynomial result;
	result.insert(ComplexMonomial::BuildH(id), value, 1);
	return result;
}
ComplexPolynomial ComplexPolynomial::BuildH(PosInt id, ComplexNumber const & value) {
	ComplexPolynomial result;
	result.insert(ComplexMonomial::BuildH(id), value, 1);
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

int ComplexPolynomial::degree()const {
	return _degree->empty() ? 0 : _degree->rbegin()->first;
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

void ComplexPolynomial::insert(ComplexMonomialPtr const & lhs, ComplexNumber rhs, ComplexNumber factor) {
	auto result = _terms->insert({ lhs, factor*rhs });
		
	if (!result.second) {
		result.first->second += factor*rhs;

		if (std::abs(result.first->second) > 0) {

		}
		else {
			_terms->erase(result.first);
			_degree->find(lhs->degree())->second.erase(lhs);
		}
	}
	else {
		(*_degree)[lhs->degree()].insert(lhs);
	}

}

void ComplexPolynomial::insert(ComplexTerms::value_type const &term, ComplexNumber factor) {
	insert(term.first, term.second, factor);
}

void ComplexPolynomial::insert(ComplexPolynomial const & rhs, ComplexNumber factor) {
	if (std::abs(factor) != 0) {
		for (auto const & term : rhs.terms()) {
			insert(term, factor);
		}
	}
}
std::ostream & ComplexPolynomial::print(std::ostream & stream) const {
	//for (auto const & term : terms()) {
	//	stream << format(term.second);
	//	stream << term.first;
	//}
	for (auto const & term : terms()) {
		if (term.first->degree() > 0) {
			stream << format(term.second);
			stream << term.first;
		}
		else if (term.second.imag() == 0)
			stream << term.second.real();
		else  if (term.second.real() == 0)
			stream << term.second.imag();
		else
			stream << term.second;
	}
	return stream;
}
std::ostream & ComplexPolynomial::print(std::ostream & stream, PolynomialOptimizationProblem const & rhs) const {
	for (auto const & term : terms()) {
		if (term.first->degree() > 0) {
			stream << format(term.second);
			term.first->print(stream, rhs);
		}
		else if (term.second.imag() == 0)
			stream << term.second.real();
		else  if (term.second.real() == 0)
			stream << term.second.imag();
		else
			stream << term.second;

	}
	return stream;
}

void ComplexPolynomial::operator+=(ComplexPolynomial const & rhs) {
	//std::cout << *this << " += " << rhs ;
	ComplexPolynomial result = (*this) + rhs;
	_terms = result._terms;
	_degree = result._degree;
	//std::cout << " = " << *this << std::endl;
}

void ComplexPolynomial::operator-=(ComplexPolynomial const & rhs) {
	ComplexPolynomial result = (*this) - rhs;
	_terms = result._terms;
	_degree = result._degree;
}

void ComplexPolynomial::operator*=(ComplexPolynomial const & rhs) {
	ComplexPolynomial result = (*this) * rhs;
	_terms = result._terms;
	_degree = result._degree;
}

void ComplexPolynomial::operator/=(ComplexPolynomial const & rhs) {
	ComplexPolynomial result = (*this) / rhs;
	_terms = result._terms;
	_degree = result._degree;
}

