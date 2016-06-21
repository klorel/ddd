#pragma once

#include <map>
#include <set>
#include <list>
#include <vector>
#include <map>
#include <string>

#include <iostream>
#include <iomanip>
#include <ostream>
#include <sstream>

#include <cmath>
#include <memory>
#include <algorithm>
#include <random>
#include <complex>
#include <stdexcept>
#include <functional>

typedef double Number;

typedef std::vector<Number> NumberVector;
typedef std::vector<size_t> PosIntVector;

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



typedef std::vector<int, int> Pattern;


typedef std::set<int> IntSet;
typedef std::shared_ptr<IntSet> IntSetPtr;

typedef std::vector<IntSet> SparsityPattern;

typedef std::multimap<int, int, std::greater<int> > Labels;
typedef Labels::iterator ItLabel;
typedef std::vector<ItLabel> ItLabels;



#include "algebra.h"

