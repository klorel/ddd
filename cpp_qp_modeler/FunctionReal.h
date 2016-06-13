#pragma once

#include "common.h"

class FunctionReal{
	friend class Problem;
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


	void print(std::ostream & stream, Problem const &)const;
	void print(std::ostream & stream)const;
private:
	mutable NumberPtr _c;
	mutable LinearTermPtr _l;
	mutable QuadraticTermPtr _q;
private:
	void operator+=(FunctionReal const &);
	void operator-=(FunctionReal const &);
	void operator*=(FunctionReal const &);
	void operator/=(FunctionReal const &);
private:
	void allocate();

	void add(Number value, Number factor);
	void add(size_t index, Number value, Number factor);
	void add(size_t index1, size_t index2, Number value, Number factor);

	void add(LinearTerm::value_type const &, Number factor);
	void add(QuadraticTerm::value_type const &, Number factor);
	void add(Index2 const &, Number value, Number factor);
};

