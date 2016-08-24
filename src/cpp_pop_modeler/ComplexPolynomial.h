#pragma once

#include "ComplexMonomialPredicate.h"

class ComplexPolynomial {
public:
	friend ComplexPolynomial operator+(ComplexPolynomial const & rhs);
	friend ComplexPolynomial operator-(ComplexPolynomial const & rhs);

	friend ComplexPolynomial operator+(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
	friend ComplexPolynomial operator-(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
	friend ComplexPolynomial operator*(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
	friend ComplexPolynomial operator/(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
public:
	ComplexPolynomial H()const;
	ComplexPolynomial conjugate()const;
	void get_all_monomial(ComplexMonomialPtr2Int & output)const;
	void clear();
public:
	ComplexPolynomial();
	ComplexPolynomial(Number value);
	ComplexPolynomial(Number real, Number imag);
	ComplexPolynomial(ComplexNumber const & value);
	ComplexPolynomial(ComplexMonomialPtr value);

	static ComplexPolynomial i();
	static ComplexPolynomial Build(PosInt id);
	static ComplexPolynomial Build(PosInt id, Number value);
	static ComplexPolynomial Build(PosInt id, ComplexNumber const & value);

	static ComplexPolynomial BuildH(PosInt id);
	static ComplexPolynomial BuildH(PosInt id, Number value);
	static ComplexPolynomial BuildH(PosInt id, ComplexNumber const & value);
	static std::vector<ComplexPolynomial> BuildVector(PosInt size);
	static std::vector<ComplexPolynomial> BuildVectorH(PosInt size);

public:
	void operator+=(ComplexPolynomial const & rhs);
	void operator-=(ComplexPolynomial const & rhs);
	void operator*=(ComplexPolynomial const & rhs);
	void operator/=(ComplexPolynomial const & rhs);

	ComplexTerms const & terms() const;
	int degree()const;
	bool isConstant()const;
	ComplexNumber constant()const;

	void insert(ComplexPolynomial const & rhs, ComplexNumber factor);
	void insert(ComplexMonomialPtr const & rhs, ComplexNumber, ComplexNumber factor);
	void insert(ComplexTerms::value_type const &, ComplexNumber factor);
	std::ostream & print(std::ostream & stream) const;
	std::ostream & print(std::ostream & stream, PolynomialOptimizationProblem const & rhs) const;

private:
	ComplexTermsPtr _terms;
	Degree2TermsPtr _degree;
};


