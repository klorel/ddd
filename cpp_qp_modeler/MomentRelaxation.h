#pragma once

#include "common.h"

class MomentRelaxation {
public:

public:
	MomentRelaxation(size_t nvariables, size_t order);
	~MomentRelaxation();
public:
	double size()const;
	
private:
	size_t _order;
	size_t _nvariables;
};