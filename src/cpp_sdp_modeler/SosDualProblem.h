
#pragma once

#include "common.h"
#include "PolynomialOptimizationProblem.h"
#include "SdpProblem.h"
#include "ConicProblem.h"

typedef void * MSKenv_t;
typedef void * MSKtask_t;

typedef std::map<ComplexMonomialPtr, IntPair2Dbl, ComplexMonomialPredicate> MonomialDecomposition;
class SosDualProblem {
public:
	SosDualProblem(PolynomialOptimizationProblem &);
	~SosDualProblem();
	void set_up_moment(int order);
	void run(int order);

	void build(int order);

public:
	int numcon()const;
	int numbarvar()const;
	int lenbarvar()const;
	void add_obj(MSKenv_t & env, MSKtask_t task);
	void add_ctr(MSKenv_t & env, MSKtask_t task, int id_ctr);
	void solve(MSKenv_t & env, MSKtask_t task);
	std::ostream & print(std::ostream &)const;

	void build_obj();
	void build_ctr(int id_ctr);

private:
	// B_alpha (0) and C_i_alpha (i>0)
	Matrix4 _matrix;

	PolynomialOptimizationProblem * _pop;
	std::vector<ComplexMonomialPtr> _id2monomial;
	ComplexMonomialPtr2Int _monomial2id;
	MonomialDecomposition _B_alpha;
	std::vector<MonomialDecomposition> C_i_alpha;
	NumberVector _rhs;
	NumberVector _c;
	
	int _n;
	NumberVector _y;
	IntVector _lenbarvar;
public:
	ConicProblem _conic;
};
