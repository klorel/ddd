
#pragma once

#include "common.h"
class SdpProblem;
void get_sdp_1(SdpProblem & sdp);

class SdpProblem {
public:
	typedef std::array<int, 4> Coef;
	typedef std::map<Coef, Number> Matrix;
	class Block
	{
	public:
		int _size;
		int _begin;
	public:
		int nz()const {
			return _size < 0 ? std::abs(_size) : _size*(_size + 1) / 2;
		}
	};
	typedef std::vector<Block> Blocks;
public:
	int newBlock(int);

	void add(int, int, int, int, Number);
	
	int newy(double cost);
	
	void print(std::ostream & stream)const;
	void print(std::string const &)const;

	void clear();

	int nvars()const;

	std::ostream & print(std::ostream &, std::string const & ) const;

	int id(int i, int j, int k)const;
	int nz()const;
	int nblock()const;
	int nctr()const;
public:
	NumberVector _b;
	Matrix _matrix;
	Blocks _blocks;
};

int numBefore(int n, int line);
