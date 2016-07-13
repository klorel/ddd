#pragma once
#include "common.h"

class RealMonomial;
class ComplexMonomial;
typedef std::shared_ptr<ComplexMonomial> ComplexMonomialPtr;
typedef std::shared_ptr<RealMonomial> RealMonomialPtr;
std::ostream & operator<<(std::ostream & stream, ComplexMonomial const & rhs);

#define __DEGUB_ORDERING__ 0

class RealMonomial {
public:

public:
	PosInt2PosInt & alpha() { return _alpha; }
	PosInt2PosInt const & alpha()const { return _alpha; }
	
	RealMonomial() :_alpha(), _degree(0) {

	}
	PosInt degree()const { return _degree; }
	PosInt &degree() { return _degree; }

	std::ostream & print(std::ostream & stream, bool const & isH)const {
		for (auto const & alpha : alpha()) {
			stream << "z" << (isH ? "H" : "") << "[" << alpha.first << "]";
			if (alpha.second > 1)
				stream << "^" << alpha.second;
		}
		return stream;
	}
private:
	PosInt2PosInt _alpha;
	PosInt _degree;
private:
};
inline std::ostream & operator<<(std::ostream & stream, RealMonomial const & rhs) {
	rhs.print(stream, true);
	return stream;
}
inline bool operator<(RealMonomial const &  lhs, RealMonomial const & rhs) {
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
	std::cout << "RealMonomial "<<lhs << (result ? "<" : "=>") << rhs << std::endl;
#endif
	return result;
}

inline bool operator>(RealMonomial const &  lhs, RealMonomial const & rhs) {
	bool result;
	if (lhs.degree() == rhs.degree())
		result = lhs.alpha() > rhs.alpha();
	else
		result = lhs.degree() > rhs.degree();
#if __DEGUB_ORDERING__
	std::cout << "RealMonomial " << lhs << (result ? ">" : "<=") << rhs << std::endl;
#endif
	return result;
}

inline bool operator==(RealMonomial const &  lhs, RealMonomial const & rhs) {
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

class ComplexMonomial {
	enum Z_ORDER {
		Z,
		ZH,
	};
public:
	RealMonomial & zH() { return *std::get<Z>(_non_zero); }
	RealMonomial const & zH()const { return *std::get<Z>(_non_zero); }

	RealMonomial & z() { return *std::get<ZH>(_non_zero); }
	RealMonomial const & z()const { return *std::get<ZH>(_non_zero); }

	bool zero() { return zH().degree() == 0 && z().degree() == 0; }

	std::ostream & print(std::ostream & stream)const {
		zH().print(stream, true);
		z().print(stream, false);

		return stream;
	}
	ComplexMonomialPtr conjugate()const {
		ComplexMonomialPtr result(new ComplexMonomial);
		result->_non_zero = _non_zero;
		std::swap(std::get<Z>(result->_non_zero), std::get<ZH>(result->_non_zero));
		return result;
	}
private:
	std::tuple<RealMonomialPtr, RealMonomialPtr> _non_zero;
public:
	ComplexMonomial() {
		std::get<ZH>(_non_zero) = RealMonomialPtr(new RealMonomial);
		std::get<Z>(_non_zero) = RealMonomialPtr(new RealMonomial);
	}
public:
	static ComplexMonomialPtr ZeroPtr;
	static ComplexMonomialPtr Build(PosInt);
	static ComplexMonomialPtr BuildH(PosInt);
};
inline std::ostream & operator<<(std::ostream & stream, ComplexMonomial const & rhs) {
	rhs.print(stream);
	return stream;
}
inline std::ostream & operator<<(std::ostream & stream, ComplexMonomialPtr const & rhs) {
	return stream << *rhs;
}
class ComplexMonomialPredicate {
public:
	bool operator()(ComplexMonomial const & lhs, ComplexMonomial const & rhs)const {
		bool result;
		if (lhs.zH() == rhs.zH()) {
			result = lhs.z()< rhs.z();			
		}
		else {
			result = lhs.zH() < rhs.zH();
		}
#if __DEGUB_ORDERING__
		std::cout << "ComplexMonomial " << lhs ;
		std::cout << (result ? " < " : " => ");
		std::cout << rhs << std::endl;
#endif
		return result;		
	}
	bool operator()(ComplexMonomialPtr const & lhs, ComplexMonomialPtr const & rhs) const {
		return operator()(*lhs, *rhs);
	}
};
typedef std::shared_ptr<ComplexMonomial> ComplexMonomialPtr;
typedef std::map<ComplexMonomialPtr, ComplexNumber, ComplexMonomialPredicate> ComplexTerms;
typedef std::shared_ptr<ComplexTerms> ComplexTermsPtr;

class ComplexPolynomial {
public:
	friend ComplexPolynomial operator+(ComplexPolynomial const & rhs);
	friend ComplexPolynomial operator-(ComplexPolynomial const & rhs);

	friend ComplexPolynomial operator+(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
	friend ComplexPolynomial operator-(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
	friend ComplexPolynomial operator*(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
	friend ComplexPolynomial operator/(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
public:
	ComplexPolynomial conjugate()const {
		ComplexPolynomial result;		
		for (ComplexTerms::const_reverse_iterator rit(_terms->rbegin()); rit != _terms->rend(); ++rit) {
			(*result._terms)[rit->first->conjugate()] = std::conj(rit->second);
		}
		return result;
	}
public:
	ComplexPolynomial() :_terms(new ComplexTerms) {

	}
	ComplexPolynomial(Number value) :_terms(new ComplexTerms) {
		terms()[ComplexMonomial::ZeroPtr] = value;
	}
	ComplexPolynomial(Number real, Number imag) :_terms(new ComplexTerms) {
		terms()[ComplexMonomial::ZeroPtr] = ComplexNumber(real, imag);
	}
	ComplexPolynomial(ComplexNumber const & value) :_terms(new ComplexTerms) {
		terms()[ComplexMonomial::ZeroPtr] = value;
	}
	static ComplexPolynomial i() {
		ComplexPolynomial result(0,1);
		return result;

	}
	static ComplexPolynomial Build(PosInt id) {
		ComplexPolynomial result;
		result.terms()[ComplexMonomial::Build(id)] = 1;
		return result;
	}
	static ComplexPolynomial Build(PosInt id, Number value) {
		ComplexPolynomial result;
		result.terms()[ComplexMonomial::Build(id)] = value;
		return result;
	}
	static ComplexPolynomial Build(PosInt id, ComplexNumber const & value) {
		ComplexPolynomial result;
		result.terms()[ComplexMonomial::Build(id)] = value;
		return result;
	}

	static ComplexPolynomial BuildH(PosInt id) {
		ComplexPolynomial result;
		result.terms()[ComplexMonomial::BuildH(id)] = 1;
		return result;
	}
	static ComplexPolynomial BuildH(PosInt id, Number value){
		ComplexPolynomial result;
		result.terms()[ComplexMonomial::BuildH(id)] = value;
		return result;
	}
	static ComplexPolynomial BuildH(PosInt id, ComplexNumber const & value) {
		ComplexPolynomial result;
		result.terms()[ComplexMonomial::BuildH(id)] = value;
		return result;
	}
	static std::vector<ComplexPolynomial> BuildVector(PosInt size) {
		std::vector<ComplexPolynomial> result(size);
		for (PosInt i(0); i < size; ++i)
			result[i] = Build(i);
		return result;
	}
	static std::vector<ComplexPolynomial> BuildVectorH(PosInt size) {
		std::vector<ComplexPolynomial> result(size);
		for (PosInt i(0); i < size; ++i)
			result[i] = BuildH(i);
		return result;
	}
private:
	ComplexTerms & terms() { return *_terms; }

public:
	ComplexTerms const & terms() const { return *_terms; }

	bool isConstant()const {
		return _terms->size() == 1 && _terms->begin()->first->zero();
	}
	ComplexNumber constant()const {
		return _terms->begin()->second;
	}

	void insert(ComplexPolynomial const & rhs, ComplexNumber factor) {
		for (auto const & term : rhs.terms()) {
			terms()[term.first] += factor*term.second;
		}
	}
	std::ostream & print(std::ostream & stream) const {
		for (auto const & term : terms()) {
			stream << format(term.second);
			stream << term.first;
		}
		return stream;
	}
private:
	ComplexTermsPtr _terms;
};

inline std::ostream & operator<<(std::ostream & stream, ComplexPolynomial const & rhs) {
	rhs.print(stream);
	return stream;
}

ComplexPolynomial operator+(ComplexPolynomial const & rhs);
ComplexPolynomial operator-(ComplexPolynomial const & rhs);

ComplexPolynomial operator+(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
ComplexPolynomial operator-(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
ComplexPolynomial operator*(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
ComplexPolynomial operator/(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);

inline ComplexPolynomial operator+(ComplexPolynomial const & rhs) {
	ComplexPolynomial result;
	result.insert(rhs, +1);
	return result;
}
inline ComplexPolynomial operator-(ComplexPolynomial const & rhs) {
	ComplexPolynomial result;
	result.insert(rhs, -1);
	return result;
}

inline ComplexPolynomial operator+(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs) {
	ComplexPolynomial result;
	result.insert(lhs, +1);
	result.insert(rhs, +1);
	return result;
}
inline ComplexPolynomial operator-(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs) {
	ComplexPolynomial result;
	result.insert(lhs, +1);
	result.insert(rhs, -1);
	return result;
}
inline ComplexPolynomial operator*(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs) {
	ComplexPolynomial result;
	
	return result;
}
inline ComplexPolynomial operator/(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs) {
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
			result.insert(lhs, ComplexNumber(1,0) / constant);
		}
	}
	return result;
}