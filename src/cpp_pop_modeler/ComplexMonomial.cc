#include "ComplexMonomial.h"

#include "PolynomialOptimizationProblem.h"

ComplexMonomialPtr  ComplexMonomial::ZeroPtr = ComplexMonomialPtr(new ComplexMonomial);


ComplexMonomialPtr ComplexMonomial::Build(int id) {
	ComplexMonomialPtr result(new ComplexMonomial);
	result->_id2degree[id] = 1;
	result->_degreeToid[1].insert(id);
	return result;
}
ComplexMonomialPtr ComplexMonomial::BuildH(int id) {
	ComplexMonomialPtr result(new ComplexMonomial);
	result->_id2degree[id] = -1;
	result->_degreeToid[-1].insert(id);
	return result;
}

int ComplexMonomial::degree()const {
	return zero() ? 0 : std::max(std::abs(_degreeToid.begin()->first), std::abs(_degreeToid.rbegin()->first));
}

bool ComplexMonomial::zero() const {
	return _degreeToid.empty();
}

std::ostream & ComplexMonomial::print(std::ostream & stream)const {
	if (!zero()) {
		for (auto const & alpha : _id2degree) {
			stream << "z" << "[" << alpha.first << "]" << (alpha.second < 0 ? "H" : "");
			if (std::abs(alpha.second) > 1)
				stream << "^" << std::abs(alpha.second);
		}
	}
	else {
		stream << "1";
	}
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
	for (Int2IntSet::const_reverse_iterator rit(_degreeToid.rbegin()); rit != _degreeToid.rend(); ++rit) {
		result->_degreeToid[-rit->first] = rit->second;
	}
	return result;
}

ComplexMonomial::ComplexMonomial() {
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
	for (auto const & kvp : lhs._id2degree) {
		result._id2degree[kvp.first] += kvp.second;
	}
	for (auto const & kvp : rhs._id2degree) {
		result._id2degree[kvp.first] += kvp.second;
	}
	for (Int2Int::iterator it(result._id2degree.begin()); it != result._id2degree.end(); ) {
		if (it->second == 0) {
			result._id2degree.erase(it++);
		}
		else {
			result._degreeToid[it->second].insert(it->first);
			++it;
		}
	}
	return ptr;
}