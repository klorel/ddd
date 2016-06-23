
#pragma once

#include "common.h"
class Pooling;
class Pooling {
public:
	class Structure{
	public:
		IntSet S;
		IntSet I;
		IntSet T;
		IntPairSet SI;
		IntPairSet IT;
		IntPairSet ST;
		IntPairSet TK;
		IntPairSet IK;
		IntTripleSet SIT;
		IntSet bL;
		IntSet bU;

		void build(Pooling const &);
	};
public:
	void pqFormulation(Problem & result)const;
public:
	Pooling(size_t n = 0, size_t k = 0);
	virtual ~Pooling();

	void newArc(int, int);
	void get(Structure &) const;
	
	Number c(int, int) const;
	Number q(int, int) const;

	Number c(IntPair const &) const;
	Number q(IntPair const &) const;

	Number bL(int) const;
	Number bU(int) const;

	void read(std::string const &);
	void out(std::ostream & = std::cout) const;
public:
	size_t n() const;
	size_t k() const;

	IntSet const & nPlus(int i) const;
	IntSet const & nMoins(int i) const;
	void outAmpl(std::ostream & stream = std::cout) const;
private:
	void allocate(size_t n, size_t k);
protected:
	// number of nodes
	size_t _n;
	// number of attribute
	size_t _k;
	//
	std::vector<IntSet> _nPlus;
	std::vector<IntSet> _nMoins;
	// capacity
	NumberDenseMatrix _b;
	// attibute upper bound / rate
	NumberDenseMatrix _q;
	// cost
	NumberDenseMatrix _c;
};
