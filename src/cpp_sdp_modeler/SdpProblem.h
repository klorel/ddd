
#pragma once

#include "common.h"
class SdpProblem;
void get_sdp_1(SdpProblem & sdp);


class SdpProblem {
public:
	class Block
	{
	public:
		int _size;
		int _begin;
		int _previous_dim;
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
	int dim()const;
	int nblock()const;
	int nctr()const;
	void dual(Matrix4 &)const;

	void read(std::string const &);
	void matrix_completion(IntSetPtrSet & output)const;
	void sparsity_pattern_1(SparsityPattern & output)const;
	void sparsity_pattern_2(SparsityPattern & output)const;
public:
	NumberVector _b;
	Matrix4 _matrix;
	Blocks _blocks;
};

int numBefore(int n, int line);
