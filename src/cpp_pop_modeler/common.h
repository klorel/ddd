#pragma once

#include <map>
#include <set>
#include <list>
#include <vector>
#include <map>
#include <string>
#include <array>

#include <iostream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <fstream>

#include <cmath>
#include <memory>
#include <algorithm>
#include <random>
#include <complex>
#include <stdexcept>
#include <functional>

typedef double Number;

typedef std::vector<Number> NumberVector;

typedef std::vector<NumberVector> NumberDenseMatrix;

typedef std::vector<size_t> PosIntVector;
typedef std::vector<int> IntVector;

typedef std::complex<double> ComplexNumber;
typedef std::vector<ComplexNumber> ComplexNumbers;
typedef std::shared_ptr<ComplexNumber> ComplexNumberPtr;
typedef std::pair<int, int> Index2;

typedef std::map<int, Number> LinearTerm;
typedef std::map<Index2, Number> QuadraticTerm;

typedef std::vector<std::string> StrVector;
typedef std::shared_ptr<std::string> StrPtr;
typedef std::vector<StrPtr> StrPtrVector;
typedef std::map<std::string, int> Str2Int;

typedef std::shared_ptr<int>  IntPtr;
typedef std::shared_ptr<Number>  NumberPtr;
typedef std::shared_ptr<LinearTerm> LinearTermPtr;
typedef std::shared_ptr<QuadraticTerm> QuadraticTermPtr;

Index2 get_index(int i, int j);

bool isZero(Number value);
bool isZero(NumberPtr const & value);
bool isZero(ComplexNumber value);

std::string format(Number value);
std::string format(ComplexNumber const & value);

class Problem;
class FunctionReal;
class FunctionComplex;
class Constraint;

typedef std::vector<Constraint> Constraints;


typedef unsigned int PosInt;
typedef std::map<PosInt, PosInt> PosInt2PosInt;
typedef std::shared_ptr<PosInt2PosInt> PosInt2PosIntPtr;
typedef std::shared_ptr<PosInt> PosIntPtr;

typedef std::vector<int, int> Pattern;
typedef std::map<int, int> Int2Int;
typedef std::shared_ptr<Int2Int> Int2IntPtr;
typedef std::vector<Int2IntPtr> Alphas;
typedef std::pair<int, int> IntPair;
typedef std::pair<IntPair, int> IntTriple;
typedef std::map<IntPair, int> IntPair2Int;
typedef std::set<IntPair> IntPairSet;
typedef std::set<IntTriple> IntTripleSet;


typedef std::set<int> IntSet;
typedef std::shared_ptr<IntSet> IntSetPtr;

typedef std::vector<IntSet> SparsityPattern;

typedef std::multimap<int, int, std::greater<int> > Labels;
typedef Labels::iterator ItLabel;
typedef std::vector<ItLabel> ItLabels;

class IntSetPredicate {
public:
	bool operator()(IntSetPtr const & lhs, IntSetPtr const & rhs)const {
		return *lhs < *rhs;
	}
};

typedef std::set<IntSetPtr, IntSetPredicate> IntSetPtrSet;

enum Sense {
	LEQ, EQ, GEQ, RNG
};

Number posInfinity();
Number negInfinity();


#include <SparseCholesky>

typedef Eigen::SparseMatrix<double, 0, int> SparseMatrix;
typedef Eigen::Triplet<double> Triplet;
typedef std::vector<Triplet> Triplets;

// forward declarations
class RealMonomial;
class ComplexMonomial;
class ComplexPolynomial;
class ComplexMonomialPredicate;

// typedef
typedef std::shared_ptr<ComplexMonomial> ComplexMonomialPtr;
typedef std::shared_ptr<RealMonomial> RealMonomialPtr;


typedef std::shared_ptr<ComplexMonomial> ComplexMonomialPtr;
typedef std::map<ComplexMonomialPtr, ComplexNumber, ComplexMonomialPredicate> ComplexTerms;
typedef std::shared_ptr<ComplexTerms> ComplexTermsPtr;

// operator overloading
std::ostream & printAlpha(std::ostream & stream, IntVector const & rhs);
std::ostream & operator<<(std::ostream &, Problem const &);
std::ostream & operator<<(std::ostream &, Constraint const &);

std::ostream & operator<<(std::ostream & stream, RealMonomial const & rhs);
std::ostream & operator<<(std::ostream & stream, ComplexMonomial const & rhs);
std::ostream & operator<<(std::ostream & stream, ComplexMonomialPtr const & rhs);
std::ostream & operator<<(std::ostream & stream, ComplexPolynomial const & rhs);

RealMonomialPtr operator+(RealMonomial const & lhs, RealMonomial const & rhs);

bool operator<(RealMonomial const &  lhs, RealMonomial const & rhs);
bool operator>(RealMonomial const &  lhs, RealMonomial const & rhs);
bool operator==(RealMonomial const &  lhs, RealMonomial const & rhs);

ComplexMonomialPtr operator+(ComplexMonomial const & lhs, ComplexMonomial const & rhs);


ComplexPolynomial operator+(ComplexPolynomial const & rhs);
ComplexPolynomial operator-(ComplexPolynomial const & rhs);

ComplexPolynomial operator+(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
ComplexPolynomial operator-(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
ComplexPolynomial operator*(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);
ComplexPolynomial operator/(ComplexPolynomial const & lhs, ComplexPolynomial const & rhs);

void operator<<(Problem & problem, ComplexPolynomial const & rhs);

Constraint operator<=(ComplexPolynomial const & lhs, ComplexNumber const & rhs);
Constraint operator==(ComplexPolynomial const & lhs, ComplexNumber const & rhs);
Constraint operator>=(ComplexPolynomial const & lhs, ComplexNumber const & rhs);


Constraint operator<=(ComplexNumber const & lhs, ComplexPolynomial const & rhs);
Constraint operator==(ComplexNumber const & lhs, ComplexPolynomial const & rhs);
Constraint operator>=(ComplexNumber const & lhs, ComplexPolynomial const & rhs);

template<class T1> std::string Str(T1 const & t1) {
	std::stringstream buffer;
	buffer << t1;
	return buffer.str();
}
template<class T1, class T2> std::string Str(T1 const & t1, T2 const & t2) {
	std::stringstream buffer;
	buffer << t1 << t2;
	return buffer.str();
}
template<class T1, class T2, class T3> std::string Str(T1 const & t1, T2 const & t2, T3 const & t3) {
	std::stringstream buffer;
	buffer << t1 << t2 << t3;
	return buffer.str();
}
template<class T1, class T2, class T3, class T4> std::string Str(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4) {
	std::stringstream buffer;
	buffer << t1 << t2 << t3<<t4;
	return buffer.str();
}

#define __DEGUB_ORDERING__ 0

