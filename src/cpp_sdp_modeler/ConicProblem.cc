#include "ConicProblem.h"

#include "mosek.h"


IndexedPool2Square & ConicProblem::new_block(int  n, bool is_sdp) {
	std::string name = Str(is_sdp ? "SDP_X" : "LP_X_", is_sdp ? _sdp_blocks.size() : _lp_blocks.size());
	return new_block(name, n, is_sdp);
}

IndexedPool2Square & ConicProblem::new_block(std::string const & name, int  n, bool is_sdp) {
	IndexedPool2SquarePtr ptr(new IndexedPool2Square(name, is_sdp ? _sdp_names : _lp_names, n, is_sdp ?(int) _sdp_blocks.size() : (int)_lp_blocks.size()));
	if (is_sdp)
		_sdp_blocks.push_back(ptr);
	else
		_lp_blocks.push_back(ptr);
	return *ptr;
}

int ConicProblem::new_ctr(Number rhs) {
	_rhs.push_back(rhs);
	return (int) _rhs.size();
}
Number & ConicProblem::rhs(int idctr) {
	return _rhs[idctr - 1];
}
Number ConicProblem::rhs(int idctr) const {
	return _rhs[idctr - 1];
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
	build_ctr(env, task);
	build_lp_var(env, task);
	build_sdp_var(env, task);
	build_lp_mat(env, task);
	build_sdp_mat(env, task);

}
void ConicProblem::build_ctr(MSKenv_t & env, MSKtask_t & task) {
	std::vector<MSKboundkeye> fx(_rhs.size(), MSK_BK_FX);
	IntVector indexes(fx.size());
	for (int i(0); i < fx.size(); ++i)
		indexes[i] = i;

	MSKrescodee  result;
	result = MSK_appendcons(task, (int)_rhs.size());
	if (result != MSK_RES_OK)throw std::invalid_argument("MSK_appendcons ");

	result = MSK_putconboundlist(task, (int)fx.size(), indexes.data(), fx.data(), _rhs.data(), _rhs.data());
	if (result != MSK_RES_OK)throw std::invalid_argument("MSK_putconboundlist ");
}


void ConicProblem::build_sdp_var(MSKenv_t & env, MSKtask_t & task) {
	if (!_sdp_blocks.empty()) {
		IntVector dimbarvar;
		for (int i(0); i < _sdp_blocks.size(); ++i) {
			dimbarvar.push_back(_sdp_blocks[i]->dim());
		}
		MSKrescodee result = MSK_appendbarvars(task, (int)_sdp_blocks.size(), dimbarvar.data());
		if (result != MSK_RES_OK)throw std::invalid_argument("MSK_appendbarvars ");
	}
}
void ConicProblem::build_lp_var(MSKenv_t & env, MSKtask_t & task) {
	if (!_lp_blocks.empty()) {
		int last_len(_lp_blocks.back()->dim());
		int num = _lp_blocks.back()->id(last_len - 1, last_len - 1);
		MSKrescodee  result = MSK_appendvars(task, num + 1);
	}
}
void ConicProblem::build_lp_mat(MSKenv_t & env, MSKtask_t & task) {
	if (!_lp.empty()) {
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
}
void ConicProblem::build_sdp_mat(MSKenv_t & env, MSKtask_t & task) {
	if (!_sdp.empty()) {
		IntVector idctr(_sdp.size());
		IntVector idmat(_sdp.size());
		IntVector idfirst(_sdp.size());
		IntVector idsecond(_sdp.size());
		NumberVector val(_sdp.size());

		int n(0);
		int obj_end(0);
		for (auto const & kvp : _sdp) {
			idctr[n] = kvp.first[0] - 1;
			idmat[n] = kvp.first[1];
			idfirst[n] = kvp.first[2];
			idsecond[n] = kvp.first[3];
			val[n] = kvp.second;
			std::cout <<std::setw(6) << idctr[n];
			std::cout <<std::setw(6) << idmat[n];
			std::cout <<std::setw(6) << idfirst[n];
			std::cout << std::setw(6) << idsecond[n];
			std::cout << std::setw(15) << val[n];
			std::cout << std::endl;
			if (kvp.first[0] > 0 && obj_end == 0) {
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
}



void ConicProblem::solve(MSKenv_t & env, MSKtask_t task) {
	MSKrescodee  r = MSK_RES_OK;
	MSKrescodee trmcode;

	/* Run optimizer */
	r = MSK_optimizetrm(task, &trmcode);

	/* Print a summary containing information
	about the solution for debugging purposes*/
	MSK_solutionsummary(task, MSK_STREAM_MSG);

	//NumberVector barx;

	if (r == MSK_RES_OK)
	{
		MSKsolstae solsta;
		double primalobj;
		double dualobj;
		MSK_getsolsta(task, MSK_SOL_ITR, &solsta);

		switch (solsta)
		{
		case MSK_SOL_STA_OPTIMAL:
		case MSK_SOL_STA_NEAR_OPTIMAL:
			//for (int i(0); i < NUMBARVAR; ++i) {
			//	barx.assign(LENBARVAR[i], -1);
			//	MSK_getbarxj(task,
			//		MSK_SOL_ITR,    /* Request the interior solution. */
			//		i,
			//		barx.data());
			//for (int j(0); j < LENBARVAR[i]; ++j) {
			//	std::cout << std::setw(5) << i;
			//	std::cout << std::setw(5) << j;
			//	std::cout << std::setw(15) << barx[j] << std::endl;;
			//}
			//}
			MSK_getprimalobj(task, MSK_SOL_ITR, &primalobj);
			MSK_getprimalobj(task, MSK_SOL_ITR, &dualobj);
			std::cout << "primalobj  = " << primalobj << std::endl;
			std::cout << "dualobj    = " << dualobj << std::endl;
			//_y.assign(numcon(), 0);
			//MSK_gety(task, MSK_SOL_ITR, _y.data());
			break;
		case MSK_SOL_STA_DUAL_INFEAS_CER:
		case MSK_SOL_STA_PRIM_INFEAS_CER:
		case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
		case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
			printf("Primal or dual infeasibility certificate found.\n");
			break;

		case MSK_SOL_STA_UNKNOWN:
			printf("The status of the solution could not be determined.\n");
			break;
		default:
			printf("Other solution status.");
			break;
		}
	}
	else
	{
		printf("Error while optimizing.\n");
	}
	//for (auto const & alpha : _monomial2id) {
	//	if (alpha.second>0)
	//		std::cout << std::setw(25) << alpha.first << " = " << _y[alpha.second - 1] << std::endl;
	//}

	if (r != MSK_RES_OK)
	{
		/* In case of an error print error code and description. */
		char symname[MSK_MAX_STR_LEN];
		char desc[MSK_MAX_STR_LEN];

		printf("An error occurred while optimizing.\n");
		MSK_getcodedesc(r,
			symname,
			desc);
		printf("Error %s - '%s'\n", symname, desc);
	}
}

void ConicProblem::get_solution(MSKenv_t & env, MSKtask_t & task, ConicSolution & result) {
	
	MSKrescodee  r = MSK_RES_OK;
	MSKrescodee trmcode = MSK_RES_OK;
	
	// Print a summary containing information about the solution for debugging purposes
	MSK_solutionsummary(task, MSK_STREAM_MSG);

	if (r == MSK_RES_OK)
	{
		MSKsolstae solsta = MSK_SOL_STA_OPTIMAL;
		MSK_getsolsta(task, MSK_SOL_ITR, &solsta);

		NumberVector barx;
		int n(0);
		switch (solsta)
		{
		case MSK_SOL_STA_OPTIMAL:
		case MSK_SOL_STA_NEAR_OPTIMAL:
			result._barx.clear();
			result._bars.clear();
			result._x.clear();
			result._y.assign(_rhs.size(), 0);
			for(auto const & block :_sdp_blocks) {
				barx.assign(block->size(), 0);
				MSKint32t barvar(block->id());
				MSKint32t dimbarvar(-1);
				MSKint64t lenbarvar(-1);
				MSK_getdimbarvarj(task, barvar, &dimbarvar);
				MSK_getlenbarvarj(task, barvar, &lenbarvar);
				if (lenbarvar != block->size())throw std::invalid_argument("lenbarvar");
				if (dimbarvar != block->dim())throw std::invalid_argument("dimbarvar");
				//std::cout << std::setw(15) << "block->size() = " << block->size() << std::endl;
				//std::cout << std::setw(15) << "block->len() = " << block->len() << std::endl;
				//std::cout << std::setw(15) << "block->id() = " << block->id() << std::endl;
				//std::cout << std::setw(15) << "barx.size() = " << barx.size() << std::endl;
				//std::cout << std::setw(15) << "dimbarvar = " << dimbarvar << std::endl;
				//std::cout << std::setw(15) << "lenbarvar = " << lenbarvar << std::endl;

				MSK_getbarxj(task, MSK_SOL_ITR, block->id(), barx.data());
				n = 0;
				for (int i(0); i < block->dim(); ++i) {
					for (int j(0); j <= i; ++j, ++n) {
						if (!isZero(barx[n])) {
							result._barx[{0, block->id(), i, j}] = barx[n];
						}
					}
				}
				MSK_getbarsj(task, MSK_SOL_ITR, block->id(), barx.data());
				n = 0;
				for (int i(0); i < block->dim(); ++i) {
					for (int j(0); j <= i; ++j, ++n) {
						if (!isZero(barx[n])) {
							result._bars[{0, block->id(), i, j}] = barx[n];
						}
					}
				}
			}
			MSK_getprimalobj(task, MSK_SOL_ITR, &result._primal);
			MSK_getdualobj(task, MSK_SOL_ITR, &result._dual);
			MSK_gety(task, MSK_SOL_ITR, result._y.data());

			break;
		case MSK_SOL_STA_DUAL_INFEAS_CER:
		case MSK_SOL_STA_PRIM_INFEAS_CER:
		case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
		case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
			printf("Primal or dual infeasibility certificate found.\n");
			break;

		case MSK_SOL_STA_UNKNOWN:
			printf("The status of the solution could not be determined.\n");
			break;
		default:
			printf("Other solution status.");
			break;
		}
	}
	else
	{
		printf("Error while optimizing.\n");
	}
	//for (auto const & alpha : _monomial2id) {
	//	if (alpha.second>0)
	//		std::cout << std::setw(25) << alpha.first << " = " << _y[alpha.second - 1] << std::endl;
	//}

	if (r != MSK_RES_OK)
	{
		/* In case of an error print error code and description. */
		char symname[MSK_MAX_STR_LEN];
		char desc[MSK_MAX_STR_LEN];

		printf("An error occurred while optimizing.\n");
		MSK_getcodedesc(r,symname,desc);
		printf("Error %s - '%s'\n", symname, desc);
	}
}