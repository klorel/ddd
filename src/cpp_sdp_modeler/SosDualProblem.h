
#pragma once

#include "common.h"
#include "PolynomialOptimizationProblem.h"
#include "SdpProblem.h"


typedef void * MSKenv_t;
typedef void * MSKtask_t;


class SosDualProblem {
public:
	SosDualProblem(PolynomialOptimizationProblem &);
	~SosDualProblem();
	void set_up_moment(int order);
	void run(int order);

public:
	int numcon()const;
	int numbarvar()const;
	int lenbarvar()const;
	void add_obj(MSKenv_t & env, MSKtask_t task);
	void add_ctr(MSKenv_t & env, MSKtask_t task, int id_ctr);
	void solve(MSKenv_t & env, MSKtask_t task);
private:
	// B_alpha (0) and C_i_alpha (i>0)
	Matrix _matrix;

	PolynomialOptimizationProblem * _pop;
	std::vector<ComplexMonomialPtr> _id2monomial;
	ComplexMonomialPtr2Int _monomial2id;
	std::map<ComplexMonomialPtr, IntPairSet, ComplexMonomialPredicate> _B_alpha;
	int _n;

};
