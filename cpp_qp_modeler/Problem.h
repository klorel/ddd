#pragma once

#include "common.h"
#include "FunctionReal.h"
#include "FunctionComplex.h"
#include "Constraint.h"
class Problem
{
public:
	Problem();
	~Problem();
public:
	std::string & name(size_t);
	std::string const & name(size_t)const;

	FunctionReal variable(size_t i)const;
	FunctionReal variable(std::string const &name, size_t i)const;
	
	FunctionComplex ivariable(std::string const &name, size_t i)const;
	FunctionComplex variable(std::string const &name1, std::string const &name2, size_t i1)const;
	FunctionComplex variable(std::string const &name1, size_t i1, std::string const &name2, size_t i2)const;
	FunctionComplex variable(std::string const &name1, size_t i1, size_t i2)const;

	void variablePool(std::string const &, size_t size);
	void ivariablePool(std::string const &, size_t size);

	size_t id(std::string const & name, size_t i)const;
	size_t first(std::string const & name)const;
	void add(Constraint const &);

	FunctionReal & minimize();
	FunctionReal const & minimize()const;

	void print(std::ostream &)const;
public:
	void addSparsityPattern(SparsityPattern & sparsityPattern)const;
private:
	StrVector _names;
	Str2Int _pools;
	Constraints _constraints;
	FunctionReal _minimize;
};

void operator<<(Problem & problem, Constraint const & rhs);