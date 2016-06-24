#pragma once

#include "common.h"
#include "FunctionReal.h"
#include "FunctionComplex.h"
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

	Constraint ctr(int key);
	Constraint ctr(std::string const &name, int i);
	Constraint ctr(std::string const &name, int i1, int i2);
	Constraint ctr(std::string const &name1, IntPair const &);

	FunctionReal variable(int i)const;
	FunctionReal variable(std::string const &name, int i)const;
	FunctionReal variable(std::string const &name1, int i1, int i2)const;
	FunctionReal variable(std::string const &name1, IntPair const &)const;

	//FunctionComplex ivariable(std::string const &name, int i)const;
	//FunctionComplex variable(std::string const &name1, std::string const &name2, int i1)const;
	//FunctionComplex variable(std::string const &name1, int i1, std::string const &name2, int i2)const;
	

	void newvarpool(std::string const &, size_t size);
	void newvarpool(std::string const & poolname, IntSet const & ids);
	void newvarpool(std::string const & poolname, IntPairSet const & ids);
	
	void newctrpool(std::string const &, size_t size);
	void newctrpool(std::string const & poolname, IntSet const & ids);
	void newctrpool(std::string const & poolname, IntPairSet const & ids);


	//void ivariablePool(std::string const &, size_t size);

	int idvar(std::string const & name, int i1)const;
	int idvar(std::string const & name, int i1, int i2)const;

	int idctr(std::string const & name, int i1)const;
	int idctr(std::string const & name, int i1, int i2)const;
	
	void add(Constraint const &);

	FunctionReal & minimize();
	FunctionReal const & minimize()const;

	void print(std::ostream &)const;

	int nvars()const;
	int nctrs()const;
public:
	void addSparsityPattern(SparsityPattern & sparsityPattern)const;
	void addSupport(SparsityPattern & sparsityPattern)const;
private:
	StrPtrVector _varnames;
	StrPtrVector _ctrnames;
	Str2Pool _varpools;
	Str2Pool _ctrpools;

	Constraints _constraints;
	FunctionReal _minimize;
};

int id(Str2Pool const & pool, std::string const & name, int i1);
int id(Str2Pool const & pool, std::string const & name, int i1, int i2);

void newpool(Str2Pool & pool, StrPtrVector & names, std::string const & poolname, size_t ids);
void newpool(Str2Pool & pool, StrPtrVector & names, std::string const & poolname, IntSet const & ids);
void newpool(Str2Pool & pool, StrPtrVector & names, std::string const & poolname, IntPairSet const & ids);
