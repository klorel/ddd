#pragma once

#include <map>
#include <set>
#include <list>
#include <vector>
#include <string>

#include <iostream>
#include <ostream>
#include <sstream>



#include <algorithm>
#include <memory>
#include <complex>
#include <stdexcept>

typedef double Number;
typedef std::complex<double> ComplexNumber;
typedef std::vector<ComplexNumber> ComplexNumbers;
typedef std::pair<size_t, size_t> Index2;

typedef std::map<size_t, Number> LinearTerm;
typedef std::map<Index2, Number> QuadraticTerm;

typedef std::vector<std::string> StrVector;
typedef std::map<std::string, size_t> Str2Int;

typedef std::shared_ptr<Number>  NumberPtr;
typedef std::shared_ptr<LinearTerm> LinearTermPtr;
typedef std::shared_ptr<QuadraticTerm> QuadraticTermPtr;

Index2 get_index(size_t i, size_t j);

bool isZero(Number value);
bool isZero(NumberPtr const & value);

std::string format(Number value);

Number PosInfinity();
Number NegInfinity();


class Problem;
class FunctionReal;
class FunctionComplex;
class Constraint;

typedef std::vector<Constraint> Constraints;

std::ostream & operator<<(std::ostream &, FunctionReal const &);
std::ostream & operator<<(std::ostream &, FunctionComplex const &);
std::ostream & operator<<(std::ostream &, Problem const &);
std::ostream & operator<<(std::ostream &, Constraint const &);


enum BranchAttributes {
	Z, YOR, YOEX, PSTOR, PSTEX, NUM
};

#include "algebra.h"

