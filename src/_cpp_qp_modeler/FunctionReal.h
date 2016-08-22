#pragma once

#include "common.h"

class FunctionReal{
	friend class PolynomialOptimizationProblem;
	friend class FunctionComplex;
	friend FunctionReal operator+(FunctionReal const &);
	friend FunctionReal operator-(FunctionReal const &);

	friend FunctionReal operator+(FunctionReal const &, FunctionReal const &);
	friend FunctionReal operator-(FunctionReal const &, FunctionReal const &);
	friend FunctionReal operator*(FunctionReal const &, FunctionReal const &);
	friend FunctionReal operator/(FunctionReal const &, FunctionReal const &);

public:
	FunctionReal();
	FunctionReal(Number);
	~FunctionReal();
public:
	FunctionReal clone()const;
public:
	bool isConstant()const;
	bool isLinear()const;
	bool isNull()const;
public:
	LinearTerm const & linear() const;
	LinearTerm & linear();

	QuadraticTerm const & quadratic() const;
	QuadraticTerm & quadratic();

	Number const & constant() const;
	Number & constant();

	void clear();

	void print(std::ostream & stream, PolynomialOptimizationProblem const &)const;
	void print(std::ostream & stream)const;
	
	void addSupport(IntSet &)const;
	void addSparsityPattern(SparsityPattern &)const;
private:
	mutable NumberPtr _c;
	mutable LinearTermPtr _l;
	mutable QuadraticTermPtr _q;
public:
	void operator+=(FunctionReal const &);
	void operator-=(FunctionReal const &);
	void operator*=(FunctionReal const &);
	void operator/=(FunctionReal const &);
private:
	void allocate();

	void add(Number value, Number factor);
	void add(int index, Number value, Number factor);
	void add(int index1, int index2, Number value, Number factor);

	void add(LinearTerm::value_type const &, Number factor);
	void add(QuadraticTerm::value_type const &, Number factor);
	void add(Index2 const &, Number value, Number factor);
};

