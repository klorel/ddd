#pragma once

#include "FunctionReal.h"

class FunctionComplex{
	friend FunctionComplex operator+(FunctionComplex const & rhs);
	friend FunctionComplex operator-(FunctionComplex const & rhs);
	friend FunctionComplex operator+(FunctionComplex const &, FunctionComplex const &);
	friend FunctionComplex operator-(FunctionComplex const &, FunctionComplex const &);
	friend FunctionComplex operator*(FunctionComplex const &, FunctionComplex const &);
	friend FunctionComplex operator/(FunctionComplex const &, FunctionComplex const &);
public:
	FunctionComplex();
	FunctionComplex(Number , Number = 0);
	~FunctionComplex();
public:
	FunctionComplex clone()const;
	FunctionComplex conjugate()const;
public:
	bool isReal()const;
	bool isImag()const;
	bool isConstant()const;
	bool isLinear()const;
	bool isNull()const;
public:
	FunctionReal const & real()const;
	FunctionReal & real();

	FunctionReal const & imag()const;
	FunctionReal & imag();

	void print(std::ostream & stream, PolynomialOptimizationProblem const &)const;

	void print(std::ostream & stream)const;

	FunctionComplex module2()const;
private:
	mutable FunctionReal _real;
	mutable FunctionReal _imag;
private:
	void operator+=(FunctionComplex const &);
	void operator-=(FunctionComplex const &);
	void operator*=(FunctionComplex const &);
	void operator/=(FunctionComplex const &);

};

