#include "FunctionReal.h"
#include "Problem.h"

FunctionReal::FunctionReal() {
	allocate();
}

FunctionReal::FunctionReal(Number rhs) {
	allocate();
	add(rhs, 1);
}

FunctionReal::~FunctionReal() {
}


void FunctionReal::allocate() {
	_c = NumberPtr(new Number(0));
	_l = LinearTermPtr(new LinearTerm);
	_q = QuadraticTermPtr(new QuadraticTerm);
}

LinearTerm const & FunctionReal::linear() const {
	return *_l;
}

LinearTerm & FunctionReal::linear() {
	return *_l;
}

QuadraticTerm const & FunctionReal::quadratic() const {
	return *_q;
}

QuadraticTerm & FunctionReal::quadratic() {
	return *_q;
}

Number const & FunctionReal::constant() const {
	return *_c;
}

Number & FunctionReal::constant() {
	return *_c;
}

FunctionReal FunctionReal::clone()const {
	FunctionReal result;
	result.constant() = constant();
	result.linear() = linear();
	result.quadratic() = quadratic();
	return result;
}

bool FunctionReal::isConstant()const {
	return isLinear() && _l->empty();
}

bool FunctionReal::isLinear()const {
	return _q->empty();
}

bool FunctionReal::isNull()const {
	return isConstant() && isZero(_c);
}

void FunctionReal::operator+=(FunctionReal const & rhs) {
	constant() += rhs.constant();
	for (auto const & kvp : rhs.linear()) {
		add(kvp.first, kvp.second, 1.0);
	}
	for (auto const & kvp : rhs.quadratic()) {
		add(kvp.first, kvp.second, 1.0);
	}
}

void FunctionReal::operator-=(FunctionReal const & rhs) {
	constant() -= rhs.constant();
	for (auto const & kvp : rhs.linear()) {
		add(kvp.first, kvp.second, -1.0);
	}
	for (auto const & kvp : rhs.quadratic()) {
		add(kvp.first, kvp.second, -1.0);
	}
}

void FunctionReal::clear() {
	constant() = 0;
	linear().clear();
	quadratic().clear();
}

void FunctionReal::operator*=(FunctionReal const & rhs) {
	if (!isLinear() && !rhs.isLinear()) {
		throw std::invalid_argument("in FunctionReal::operator*=");
	}
	else {
		if (rhs.isNull() || isNull()) {
			clear();
		}
		else {
			for (auto & lkvp : linear()) {
				for (auto & rkvp : rhs.linear()) {
					add(lkvp.first, rkvp.first, lkvp.second*rkvp.second, 1.0);
				}
			}
			if (isZero(rhs.constant())) {
				linear().clear();
			}
			else {
				for (auto & kvp : linear()) {
					kvp.second *= rhs.constant();
				}
			}
			if (!isZero(constant())) {
				for (auto const & kvp : rhs.linear()) {
					add(kvp, constant());
				}
			}
			constant() *= rhs.constant();
		}


	}
}

void FunctionReal::operator/=(FunctionReal const & rhs) {
	if (rhs.isNull() || !rhs.isConstant()) {
		throw std::invalid_argument("in FunctionReal::operator/=");
	}
	else {
		constant() /= rhs.constant();
		for (auto & kvp : linear()) {
			kvp.second /= rhs.constant();
		}
		for (auto & kvp : quadratic()) {
			kvp.second /= rhs.constant();
		}
	}
}

void FunctionReal::add(Number value, Number factor) {
	constant() += value*factor;
}

void FunctionReal::add(size_t key, Number value, Number factor) {
	LinearTerm::iterator it(linear().find(key));
	if (it == linear().end()) {
		linear()[key] = value*factor;
	}
	else {
		it->second += value*factor;
		if (isZero(it->second)) {
			linear().erase(it);
		}
	}
}

void FunctionReal::add(size_t index1, size_t index2, Number value, Number factor) {
	add(get_index(index1, index2), value, factor);
}

void FunctionReal::add(LinearTerm::value_type const & kvp, Number factor) {
	add(kvp.first, kvp.second, factor);
}

void FunctionReal::add(QuadraticTerm::value_type const & kvp, Number factor) {
	add(kvp.first.first, kvp.second, factor);
}

void FunctionReal::add(Index2 const & key, Number value, Number factor) {
	QuadraticTerm::iterator it(quadratic().find(key));
	if (it == quadratic().end()) {
		quadratic()[key] = value*factor;
	}
	else {
		it->second += value*factor;
		if (isZero(it->second)) {
			quadratic().erase(it);
		}
	}
}

void FunctionReal::print(std::ostream & stream, Problem const & problem)const {
	if (constant() != 0) {
		stream << constant();
	}
	for (auto & kvp : linear())
		stream << format(kvp.second) << problem.name(kvp.first);
	for (auto & kvp : quadratic())
		stream << format(kvp.second) << problem.name(kvp.first.first) << problem.name(kvp.first.second);
}
void FunctionReal::print(std::ostream & stream)const {
	if (constant() != 0) {
		stream << constant();
	}
	for (auto & kvp : linear())
		stream << format(kvp.second) << "x[" << kvp.first << "]";
	for (auto & kvp : quadratic())
		stream << format(kvp.second) << "x[" << kvp.first.first << "]" << "x[" << kvp.first.second << "]";
}

