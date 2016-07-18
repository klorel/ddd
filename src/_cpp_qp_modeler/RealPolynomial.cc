#include "RealPolynomial.h"
#include "Problem.h"

RealPolynomial::RealPolynomial() :_terms(new RealPolynomialTerms), _degree(new int(0)) {

}
RealPolynomial::RealPolynomial(Number const & value) : _terms(new RealPolynomialTerms), _degree(new int(0)) {
	set(value);
}
RealPolynomial::RealPolynomial(NumberPtr const &value) : _terms(new RealPolynomialTerms), _degree(new int(0)) {
	set(*value);
}

RealPolynomial::~RealPolynomial() {

}
int & RealPolynomial::degree() {
	return *_degree;
}

int const & RealPolynomial::degree()const {
	return *_degree;
}
RealPolynomialTerms & RealPolynomial::terms() {
	return *_terms;
}

RealPolynomialTerms const & RealPolynomial::terms()const {
	return *_terms;
}

void RealPolynomial::add(RealPolynomialTerm const & key, Number value) {
	add(std::make_pair(key, value), 1.0);
}

void RealPolynomial::add(RealPolynomialTerms::value_type const & kvp, Number factor) {
	auto insert_result(_terms->insert(kvp));
	if (!insert_result.second) {
		insert_result.first->second += kvp.second*factor;
		if (isZero(insert_result.first->second)) {
			_terms->erase(insert_result.first);
		}
	}
	else
		insert_result.first->second *= factor;
}

void RealPolynomial::add(RealPolynomial const & rhs, Number factor) {
	for (auto const & term : rhs.terms()) {
		add(term, factor);
	}
}

void RealPolynomial::clear() {
	terms().clear();
	degree() = 0;
}

RealPolynomial RealPolynomial::clone()const {
	RealPolynomial result;
	result.degree() = degree();
	result.terms() = terms();
	return result;
}

std::ostream &  RealPolynomial::print(std::ostream & stream, Problem const & problem)const {
	for (auto const & term : terms()) {
		stream << format(term.second);
		for (auto const & alpha : *term.first.second) {
			stream << problem.name(alpha.first);
			if (alpha.second > 1)
				stream << "^" << alpha.second;
		}
	}
	return stream;
}

std::ostream &  RealPolynomial::print(std::ostream & stream)const {
	for (auto const & term : terms()) {
		stream << format(term.second);
		stream << *term.first.second;
	}
	return stream;
}


void RealPolynomial::set(int id, Number value) {
	clear();
	Int2IntPtr key(new Int2Int);
	(*key)[id] = 1;
	terms()[{1, key}] = value;
	degree() = 1;
}
void RealPolynomial::set(int id) {
	set(id, 1.0);
}
void RealPolynomial::set(Number value) {
	clear();
	Int2IntPtr key(new Int2Int);
	terms()[{0, key}] = value;
}
Number RealPolynomial::constant()const {
	Int2IntPtr key(new Int2Int);
	RealPolynomialTerms::const_iterator it(terms().find({ 0,key }));
	if (it != terms().end())
		return it->second;
	else
		return 0;
}
RealPolynomial operator+(RealPolynomial const &rhs) {
	return rhs;
}

RealPolynomial operator-(RealPolynomial const &rhs) {
	RealPolynomial result;
	result.add(rhs, -1);
	return result;
}
RealPolynomial operator+(RealPolynomial const &lhs, RealPolynomial const &rhs) {
	RealPolynomial result;
	result.add(lhs, +1);
	result.add(rhs, +1);
	return result;

}

RealPolynomial operator-(RealPolynomial const &lhs, RealPolynomial const &rhs) {
	RealPolynomial result;
	result.add(lhs, +1);
	result.add(rhs, -1);
	return result;

}

RealPolynomial operator*(RealPolynomial const &lhs, RealPolynomial const &rhs) {
	RealPolynomial result;
	result.degree() = lhs.degree() + rhs.degree();
	for (auto const & lhs_term : lhs.terms()) {
		for (auto const & rhs_term : rhs.terms()) {
			RealPolynomialTerm term({ 0,Int2IntPtr(new Int2Int) });
			for (auto const & lhs_key : *lhs_term.first.second) {
				(*term.second)[lhs_key.first] += lhs_key.second;
				term.first = std::max(term.first, (*term.second)[lhs_key.first]);
			}
			for (auto const & rhs_key : *rhs_term.first.second) {
				(*term.second)[rhs_key.first] += rhs_key.second;
				term.first = std::max(term.first, (*term.second)[rhs_key.first]);
			}
			result.add(term, lhs_term.second*rhs_term.second);
		}
	}
	return result;

}

RealPolynomial operator/(RealPolynomial const &lhs, RealPolynomial const &rhs) {
	RealPolynomial result;
	if (rhs.degree() > 0) {
		throw std::invalid_argument("in operator/(RealPolynomial const &lhs, RealPolynomial const &rhs) rhs has a non zero degree");
	}
	else {
		Number const constant(rhs.constant());
		if (isZero(constant)) {
			throw std::invalid_argument("in operator/(RealPolynomial const &lhs, RealPolynomial const &rhs): rhs is zero");
		}
		else {
			result.add(lhs, 1 / constant);
		}
	}
	return result;
}

std::ostream & operator<<(std::ostream & stream, RealPolynomial const & rhs) {
	rhs.print(stream);
	return stream;
}
std::ostream & operator<<(std::ostream & stream, RealPolynomialTerms const & rhs) {
	for (auto const & term : rhs) {
		stream << format(term.second);
		stream << term.first.second;
	}
	return stream;
}
std::ostream & operator<<(std::ostream & stream, Int2Int const & rhs) {
	for (auto const & alpha : rhs) {
		stream << "x[" << alpha.first << "]";
		if (alpha.second > 1)
			stream << "^" << alpha.second;
	}
	return stream;
}

std::ostream & operator<<(std::ostream & stream, RealPolynomialTermsPtr const & rhs) {
	stream << *rhs;
	return stream;
}