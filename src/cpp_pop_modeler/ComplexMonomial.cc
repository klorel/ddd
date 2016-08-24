#include "ComplexMonomial.h"

#include "PolynomialOptimizationProblem.h"

ComplexMonomialPtr  ComplexMonomial::ZeroPtr = ComplexMonomialPtr(new ComplexMonomial);


ComplexMonomialPtr ComplexMonomial::Build(int id) {
	ComplexMonomialPtr result(new ComplexMonomial);
	result->_id2degree[id] = 1;
	result->_degree = 1;
	return result;
}
ComplexMonomialPtr ComplexMonomial::BuildH(int id) {
	ComplexMonomialPtr result(new ComplexMonomial);
	result->_id2degree[id] = -1;
	result->_degree = 1;
	return result;
}

int ComplexMonomial::degree()const {
	return _degree;
}

bool ComplexMonomial::zero() const {
	return _degree == 0;
}

std::ostream & ComplexMonomial::print(std::ostream & stream)const {
	std::stringstream buffer;
	if (!zero()) {
		for (auto const & alpha : _id2degree) {
			buffer << "z" << "[" << alpha.first << "]" << (alpha.second < 0 ? "H" : "");
			if (std::abs(alpha.second) > 1)
				buffer << "^" << std::abs(alpha.second);
		}
	}
	else {
		buffer << "1";
	}
	stream << buffer.str();
	return stream;
}
std::ostream & ComplexMonomial::print(std::ostream & stream, PolynomialOptimizationProblem const & rhs)const {
	for (auto const & alpha : _id2degree) {
		stream << rhs.name(alpha.first) << (alpha.second < 0 ? "H" : "");
		if (std::abs(alpha.second) > 1)
			stream << "^" << std::abs(alpha.second);
	}
	return stream;
}
ComplexMonomialPtr ComplexMonomial::conjugate()const {
	ComplexMonomialPtr result(new ComplexMonomial);
	for (auto const & kvp : _id2degree) {
		result->_id2degree[kvp.first] = -kvp.second;
	}
	result->_degree = _degree;
	return result;
}

ComplexMonomial::ComplexMonomial() :_degree(0){
}


std::ostream & operator<<(std::ostream & stream, ComplexMonomial const & rhs) {
	rhs.print(stream);
	return stream;
}

std::ostream & operator<<(std::ostream & stream, ComplexMonomialPtr const & rhs) {
	return stream << *rhs;
}

ComplexMonomialPtr operator+(ComplexMonomialPtr const & lhs, ComplexMonomialPtr const & rhs) {
	return *lhs + *rhs;
}
ComplexMonomialPtr operator+(ComplexMonomial const & lhs, ComplexMonomial const & rhs) {
	// (alpha, alphaH) + (beta, betaH) = (alpha+beta, alphaH+betaH)
	ComplexMonomialPtr ptr(new ComplexMonomial);
	Int2Int & result(ptr->_id2degree);
	for (auto const & kvp : lhs._id2degree) {
		result[kvp.first] += kvp.second;
	}
	for (auto const & kvp : rhs._id2degree) {
		result[kvp.first] += kvp.second;
	}
	for (Int2Int::iterator it(result.begin()); it != result.end(); ) {
		if (it->second == 0) {
			result.erase(it++);
		}
		else {
			ptr->_degree += it->second;
			++it;
		}
	}
	return ptr;
}


bool ComplexMonomialPredicate::operator()(ComplexMonomial const & lhs, ComplexMonomial const & rhs)const {
	bool result(false);
	if (lhs.degree() == rhs.degree())
		result = lhs._id2degree < rhs._id2degree;
	else
		result = lhs.degree() < rhs.degree();
#if __DEGUB_ORDERING__
	std::cout << "ComplexMonomial " << lhs << "(" << lhs.degree() << ")";
	std::cout << (result ? " < " : " => ");
	std::cout << rhs << "(" << rhs.degree() << ")" << std::endl;
#endif
	return result;
}