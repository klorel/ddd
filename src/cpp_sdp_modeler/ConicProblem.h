
#pragma once

#include "common.h"
#include "IndexedPool.h"

typedef void * MSKenv_t;
typedef void * MSKtask_t;

class ConicProblem {
public:
	// sdp part of the problem
	Matrix4 _sdp;
	// linear part of the problem
	Matrix2 _lp;
	// second membre
	NumberVector _rhs;
	// variables (square matrix sdp or free)
	StrPtrVector _sdp_names;
	StrPtrVector _lp_names;
	//
	std::vector<IndexedPool2SquarePtr> _lp_blocks;
	std::vector<IndexedPool2SquarePtr> _sdp_blocks;

public:
	IndexedPool2Square & new_block(std::string const & name, int n, bool is_sdp);
	IndexedPool2Square & new_block(int n, bool is_sdp);

	void add_sdp(int idctr, int idmat, int i, int j, Number v);
	void add_lp(int idctr, int idmat, int i, int j, Number v);
	int new_ctr(Number);

	void build(MSKenv_t & env, MSKtask_t & task);
	void build_ctr(MSKenv_t & env, MSKtask_t & task);
	
	void build_lp_var(MSKenv_t & env, MSKtask_t & task);
	void build_sdp_var(MSKenv_t & env, MSKtask_t & task);

	void build_lp_mat(MSKenv_t & env, MSKtask_t & task);
	void build_sdp_mat(MSKenv_t & env, MSKtask_t & task);
public:
	//IntVector _size;
	//IntVector _begin;
	//IntVector_previous_dim;
};
