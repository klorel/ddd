#include "Problem.h"
#include "Constraint.h"

Problem::Problem()
{
}


Problem::~Problem()
{
}

std::string  & Problem::name(size_t key){
	return _names[key];
}
std::string const & Problem::name(size_t key)const{
	return _names[key];
}

FunctionReal Problem::variable(size_t i)const{
	FunctionReal result;
	result.add(i, 1);
	return result;
}
FunctionReal Problem::variable(std::string const &name, size_t i)const{
	FunctionReal result;
	result.add(id(name,i), 1, 1);
	return result;
}

FunctionComplex Problem::ivariable(std::string const &name, size_t i)const {
	return variable(name + "_real", name + "_imag", i);
}

FunctionComplex Problem::variable(std::string const &name1, std::string const &name2, size_t i1)const{
	FunctionComplex result;
	result.real().add(id(name1, i1), 1, 1);
	result.imag().add(id(name2, i1), 1, 1);
	return result;
}

FunctionComplex Problem::variable(std::string const &name1, size_t i1, std::string const &name2, size_t i2)const{
	FunctionComplex result;
	result.real().add(id(name1, i1), 1, 1);
	result.imag().add(id(name2, i2), 1, 1);
	return result;
}

FunctionComplex Problem::variable(std::string const &name1, size_t i1, size_t i2)const{
	FunctionComplex result;
	result.real().add(id(name1, i1), 1, 1);
	result.imag().add(id(name1, i2), 1, 1);
	return result;
}

void Problem::variablePool(std::string const & poolname, size_t size) {
	_pools[poolname] = _names.size();
	_names.resize(_names.size() + size);
	for (size_t i(0); i < size; ++i)
	{
		std::stringstream buffer;
		buffer << poolname << "[" << i << "]";
		_names[_names.size() - size + i] = buffer.str();
	}
}
void Problem::ivariablePool(std::string const & poolname, size_t size) {
	variablePool(poolname + "_real", size);
	variablePool(poolname + "_imag", size);
}

size_t Problem::id(std::string const & name, size_t i)const{
	return first(name) + i;
}

size_t Problem::first(std::string const & name)const{
	Str2Int::const_iterator it(_pools.find(name));
	if (it == _pools.end()){
		throw std::invalid_argument("pool name was not created");
		return -1;
	}
	else{
		return it->second;
	}
}

void Problem::add(Constraint const &rhs){
	_constraints.push_back(rhs);
}


FunctionReal & Problem::minimize() {
	return _minimize;
}
FunctionReal const & Problem::minimize()const {
	return _minimize;
}

void operator<<(Problem & problem, Constraint const & rhs) {
	problem.add(rhs);
}
void operator<<(Problem & problem, FunctionReal const & rhs) {
	problem.minimize() = rhs.clone();
}

void Problem::print(std::ostream &stream)const {
	stream << "minimize ";
	minimize().print(stream, *this);
	for (auto const & ctr : _constraints) {
		ctr.print(stream, *this);
		stream<< std::endl;
	}

}