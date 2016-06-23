#pragma once

#include "common.h"

class MomentGenerator {
public:

public:
	// maximum number of variable
	// maximum order
	MomentGenerator(size_t nvariables, size_t order);
	~MomentGenerator();


	void build();
	std::ostream & print(std::ostream &,Int2Int const & alpha, size_t n) const;
private:
	size_t _order;
	size_t _nvariables;

	Alphas _alphas;
};