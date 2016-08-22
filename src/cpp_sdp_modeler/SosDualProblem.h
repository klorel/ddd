
#pragma once

#include "common.h"
#include "PolynomialOptimizationProblem.h"
#include "SdpProblem.h"

class SosDualProblem {
public:
	SosDualProblem(PolynomialOptimizationProblem &);
	~SosDualProblem();

	void run();
private:
	// B_alpha (0) and C_i_alpha (i>0)
	Matrix _matrix;

	PolynomialOptimizationProblem * _pop;

};
