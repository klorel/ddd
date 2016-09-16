#include "ConicProblem.h"

#include "mosek.h"


IndexedPool2Square & ConicProblem::new_block(int  n, bool is_sdp) {	
	std::string name = Str(is_sdp ? "SDP_X" : "LP_X_", is_sdp ? _sdp_blocks.size() : _lp_blocks.size());
	return new_block(name, n, is_sdp);
}

IndexedPool2Square & ConicProblem::new_block(std::string const & name, int  n, bool is_sdp) {
	IndexedPool2SquarePtr ptr(new IndexedPool2Square(name, is_sdp? _sdp_names:_lp_names, n, is_sdp ? _sdp_blocks.size() : _lp_blocks.size()));
	if(is_sdp)
		_sdp_blocks.push_back(ptr);
	else
		_lp_blocks.push_back(ptr);
	return *ptr;
}

int ConicProblem::new_ctr(Number rhs) {
	int const result((int)_rhs.size());
	_rhs.push_back(rhs);
	return result;
}

void ConicProblem::add_lp(int idctr, int idmat, int i, int j, Number v) {
	int const min_ij(std::min(i, j));
	int const max_ij(std::max(i, j));
	int const id(_lp_blocks[idmat]->id(min_ij, max_ij));
	// attention on double pour hors diagonale
	_lp[{idctr, id}] += (min_ij == max_ij ? 1 : 2)*v;
}
void ConicProblem::add_sdp(int idctr, int idmat, int i, int j, Number v) {
	int const min_ij(std::min(i, j));
	int const max_ij(std::max(i, j));
	_sdp[{idctr, idmat, max_ij, min_ij}] += v;
}

void ConicProblem::build(MSKenv_t & env, MSKtask_t & task) {
}
void ConicProblem::build_ctr(MSKenv_t & env, MSKtask_t & task) {
	std::vector<MSKboundkeye> fx(_rhs.size(), MSK_BK_FX);
	IntVector indexes(fx.size());
	for (int i(0); i < fx.size(); ++i)
		indexes[i] = i;
	MSKrescodee  result = MSK_putconboundlist(task, (int)fx.size(), indexes.data(), fx.data(), _rhs.data(), _rhs.data());
}


void ConicProblem::build_sdp_var(MSKenv_t & env, MSKtask_t & task) {
	if (!_sdp_blocks.empty()) {
		IntVector lenbarvar(_sdp_blocks.size());
		for (int i(0); i < _sdp_blocks.size(); ++i) {
			lenbarvar.push_back(_sdp_blocks[i]->size());
		}
		MSKrescodee result = MSK_appendbarvars(task, (int)_sdp_blocks.size(), lenbarvar.data());
		if (result != MSK_RES_OK)throw std::invalid_argument("MSK_appendbarvars ");
	}
}
void ConicProblem::build_lp_var(MSKenv_t & env, MSKtask_t & task) {
	if (!_lp_blocks.empty()) {
		int last_len(_lp_blocks.back()->len());
		int num = _lp_blocks.back()->id(last_len - 1, last_len - 1);
		MSKrescodee  result = MSK_appendvars(task, num + 1);
	}
}
void ConicProblem::build_lp_mat(MSKenv_t & env, MSKtask_t & task) {
	IntVector subi;
	IntVector subj;
	NumberVector val;

	int n(0);
	int obj_end(0);
	MSKrescodee  result = MSK_RES_OK;

	result = MSK_putclist(task, obj_end, subi.data(), val.data());
	if (result != MSK_RES_OK)throw std::invalid_argument("MSK_putclist");

	result = MSK_putaijlist(task, n - obj_end, &subi[obj_end], &subj[obj_end], &val[obj_end]);
	if (result != MSK_RES_OK)throw std::invalid_argument("MSK_putaijlist");	

}
void ConicProblem::build_sdp_mat(MSKenv_t & env, MSKtask_t & task) {	
	IntVector idctr(_sdp.size()-1);
	IntVector idmat(_sdp.size());
	IntVector idfirst(_sdp.size());
	IntVector idsecond(_sdp.size());
	NumberVector val(_sdp.size());

	int n(0);
	int obj_end(0);
	for (auto const & kvp : _sdp) {
		idctr[n] = kvp.first[0];
		idmat[n] = kvp.first[1];
		idfirst[n] = kvp.first[2];
		idsecond[n] = kvp.first[3];
		val[n] = kvp.second;

		if (kvp.first[0] > 0) {
			obj_end = n;
		}
		++n;
	}

	MSKrescodee  result(MSK_RES_OK);
	result = MSK_putbarcblocktriplet(task, obj_end, idmat.data(), idfirst.data(), idsecond.data(), val.data());
	if (result != MSK_RES_OK)throw std::invalid_argument("MSK_putbarablocktriplet");

	result = MSK_putbarablocktriplet(task, n - obj_end, &idctr[obj_end], &idmat[obj_end], &idfirst[obj_end], &idsecond[obj_end], &val[obj_end]);
	if (result != MSK_RES_OK)throw std::invalid_argument("MSK_putbarablocktriplet");

}