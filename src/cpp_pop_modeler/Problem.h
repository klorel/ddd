#pragma once

#include "common.h"
#include "Constraint.h"
#include "VariablePool.h"

class Problem
{
public:
	Problem();
	~Problem();
public:
	std::string & name(int);
	std::string const & name(int)const;
	std::string & ctrname(int);
	std::string const &ctrname(int)const;

	Constraint const &ctr(int key) const;

	Constraint &ctr(int key);
	Constraint &ctr(std::string const &name, int i);
	Constraint &ctr(std::string const &name, int i1, int i2);
	Constraint &ctr(std::string const &name1, IntPair const &);

	ComplexPolynomial variable(int i)const;
	//ComplexPolynomial variable(std::string const &name, int i)const;
	//ComplexPolynomial variable(std::string const &name1, int i1, int i2)const;
	//ComplexPolynomial variable(std::string const &name1, IntPair const &)const;

	//FunctionComplex ivariable(std::string const &name, int i)const;
	//FunctionComplex variable(std::string const &name1, std::string const &name2, int i1)const;
	//FunctionComplex variable(std::string const &name1, int i1, std::string const &name2, int i2)const;
	
	IndexedPool const & varpool(std::string const &)const;
	IndexedPool const & ctrpool(std::string const &)const;

	IndexedPool const & newvarpool(std::string const &, size_t size);
	IndexedPool const & newvarpool(std::string const & poolname, IntSet const & ids);
	IndexedPool const & newvarpool(std::string const & poolname, IntPairSet const & ids);
	
	IndexedPool const & newctrpool(std::string const &, size_t size);
	IndexedPool const & newctrpool(std::string const & poolname, IntSet const & ids);
	IndexedPool const & newctrpool(std::string const & poolname, IntPairSet const & ids);
	
	//void ivariablePool(std::string const &, size_t size);

	int idvar(std::string const & name, int i1)const;
	int idvar(std::string const & name, int i1, int i2)const;

	int idctr(std::string const & name, int i1)const;
	int idctr(std::string const & name, int i1, int i2)const;
	
	void add(Constraint const &);

	ComplexPolynomial & minimize();
	ComplexPolynomial const & minimize()const;

	void print(std::ostream &)const;

	int nvars()const;
	int nctrs()const;
public:
	//void amplExport(std::string const &)const;
	//void addSparsityPattern(SparsityPattern & sparsityPattern)const;
	//void addSupport(SparsityPattern & sparsityPattern)const;
	//void removeInequality();
private:
	StrPtrVector _varnames;
	StrPtrVector _ctrnames;
	Str2Pool _varpools;
	Str2Pool _ctrpools;

	Constraints _constraints;
	ComplexPolynomial _minimize;
};

int id(Str2Pool const & pool, std::string const & name, int i1);
int id(Str2Pool const & pool, std::string const & name, int i1, int i2);

IndexedPoolPtr newpool(Str2Pool & pool, StrPtrVector & names, std::string const & poolname, size_t ids);
IndexedPoolPtr newpool(Str2Pool & pool, StrPtrVector & names, std::string const & poolname, IntSet const & ids);
IndexedPoolPtr newpool(Str2Pool & pool, StrPtrVector & names, std::string const & poolname, IntPairSet const & ids);
