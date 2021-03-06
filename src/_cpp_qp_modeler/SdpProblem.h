
#pragma once

#include "common.h"

class PolynomialOptimizationProblem;

class SdpProblem {
public:
	typedef std::array<int, 4> Coef;
	typedef std::map<Coef, Number> Matrix;
	class Block
	{
	public:
		int _size;
		int _begin;
	};
	typedef std::vector<Block> Blocks;
public:
	int newBlock(int);

	void add(int, int, int, int, Number);
	
	int newy(double cost);
	
	void print(std::ostream & stream)const;
	void print(std::string const &)const;

	void clear();

	void sdprelaxation(PolynomialOptimizationProblem const &);
	void sparsesdp(SdpProblem &);

	void addSparsityPattern(SparsityPattern & result)const;

	int nvars()const;

	std::ostream & print(std::ostream &, std::string const & ) const;
public:
	NumberVector _b;
	Matrix _matrix;
	Blocks _blocks;
};