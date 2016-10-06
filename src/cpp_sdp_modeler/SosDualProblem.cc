#include "SosDualProblem.h"
#include "MomentGenerator.h"

#include "mosek.h"



static void MSKAPI printstr_sos_dual_problem(void *handle, MSKCONST char str[])
{
	std::cout << str;
}

SosDualProblem::SosDualProblem(PolynomialOptimizationProblem & rhs) :_pop(&rhs) {

}
SosDualProblem::~SosDualProblem() {

}

int SosDualProblem::numcon()const {
	return static_cast<int>(_id2monomial.size() - 1);
}
int SosDualProblem::numbarvar()const {
	return static_cast<int>(_pop->nctrs() + 1);
}

int SosDualProblem::lenbarvar()const {
	return _n*(_n + 1) / 2;
}

void SosDualProblem::set_up_moment(int order) {
	MomentGenerator momentGenerator(_pop->nvars(), 2 * order);
	ComplexMonomialPtrList momentVector;
	momentGenerator.build(momentVector);
	for (auto const & ptr : momentVector) {
		_monomial2id[ptr] = 0;
	}
	_id2monomial.reserve(_monomial2id.size());
	for (auto & kvp : _monomial2id) {
		kvp.second = (int)_id2monomial.size();
		_id2monomial.push_back(kvp.first);
	}
	_n = 0;
	for (ComplexMonomialPtr2Int::const_iterator it(_monomial2id.begin()); it != _monomial2id.end() && it->first->degree() <= order; ++it) {
		++_n;
		for (ComplexMonomialPtr2Int::const_iterator jt(it); jt != _monomial2id.end() && jt->first->degree() <= order; ++jt) {
			_B_alpha[it->first + jt->first].insert({ { jt->second, it->second }, 1 });
		}
	}
	//std::cout << "momentVector size is " << momentVector.size() << std::endl;
	//std::cout << "n = " << _n << std::endl;
	//for (auto & kvp : _monomial2id)
	//std::cout << "y[" << kvp.second << " ] = " << kvp.first << std::endl;
	//for (auto const & b_alpha : _B_alpha) {
	//	std::cout << b_alpha.first << " : ";
	//	for (auto const & ij : b_alpha.second) {
	//		std::cout << "(" << ij.first << ", " << ij.second << ")";
	//	}
	//	std::cout << std::endl;
	//}
	std::cout << "Number of monomials     is " << _monomial2id.size() << std::endl;
	std::cout << "Highest monomial degree is " << _monomial2id.rbegin()->first->degree() << std::endl;
}
void SosDualProblem::build(int order) {
	set_up_moment(order);
	build_obj();
	// part related to the contraints
	C_i_alpha.resize(_pop->nctrs());
	for (int i(0); i < _pop->nctrs(); ++i) {
		build_ctr(i);
	}
}
void SosDualProblem::run(int order) {
	set_up_moment(order);

	MSKenv_t     env = NULL;
	MSKtask_t    task = NULL;
	//
	MSKrescodee  r = MSK_RES_OK;
	/* Create the mosek environment. */
	r = MSK_makeenv(&env, NULL);
	if (r != MSK_RES_OK)throw std::invalid_argument("MSK_makeenv ");

	/* Create the optimization task. */
	r = MSK_maketask(env, numcon(), 0, &task);
	if (r != MSK_RES_OK)throw std::invalid_argument("MSK_maketask ");

	MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr_sos_dual_problem);
	if (r != MSK_RES_OK)throw std::invalid_argument("MSK_linkfunctotaskstream ");

	//r = MSK_putintparam(task, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_NONE);
	//if (r != MSK_RES_OK)throw std::invalid_argument("MSK_IPAR_INTPNT_SCALING ");
	//double const zero(1e-6);
	//r = MSK_putdouparam(task, MSK_DPAR_INTPNT_TOL_PFEAS, zero);
	//if (r != MSK_RES_OK)throw std::invalid_argument("MSK_DPAR_INTPNT_TOL_PFEAS ");
	//r = MSK_putdouparam(task, MSK_DPAR_INTPNT_TOL_DFEAS, zero);
	//if (r != MSK_RES_OK)throw std::invalid_argument("MSK_DPAR_INTPNT_TOL_DFEAS ");
	//r = MSK_putdouparam(task, MSK_DPAR_INTPNT_TOL_REL_GAP, zero);
	//if (r != MSK_RES_OK)throw std::invalid_argument("MSK_DPAR_INTPNT_TOL_REL_GAP ");

	r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);
	if (r != MSK_RES_OK)throw std::invalid_argument("MSK_putobjsense ");
	r = MSK_appendcons(task, numcon());
	if (r != MSK_RES_OK)throw std::invalid_argument("MSK_appendcons ");

	// part related to the objective function
	add_obj(env, task);
	// part related to the contraints
	C_i_alpha.resize(_pop->nctrs());
	for (int i(0); i < _pop->nctrs(); ++i) {
		add_ctr(env, task, i);
	}
	//print(std::cout);

	solve(env, task);

	/* Delete the task and the associated data. */
	MSK_deletetask(&task);
	/* Delete the environment and the associated data. */
	MSK_deleteenv(&env);

}

void SosDualProblem::solve(MSKenv_t & env, MSKtask_t task) {
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
			_y.assign(numcon(), 0);
			MSK_gety(task, MSK_SOL_ITR, _y.data());
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
	for (auto const & alpha : _monomial2id) {
		if(alpha.second>0)
		std::cout << std::setw(25)<<alpha.first << " = " << _y[alpha.second - 1] << std::endl;
	}

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


// max -X(1,1)
// p_alpha = <X, B_alpha>
// X sdp
void SosDualProblem::add_obj(MSKenv_t & env, MSKtask_t task) {
	MSKrescodee  r = MSK_RES_OK;
	int id_barvar;
	r = MSK_getnumbarvar(task, &id_barvar);
	if (r != MSK_RES_OK)throw std::invalid_argument("MSK_getnumbarvar ");

	_lenbarvar.push_back(lenbarvar());
	r = MSK_appendbarvars(task, 1, &_lenbarvar.back());
	if (r != MSK_RES_OK)throw std::invalid_argument("MSK_appendbarvars ");
	IntVector i;
	IntVector j;
	IntVector k;
	IntVector l;
	NumberVector v;

	j.push_back(id_barvar);
	k.push_back(0);
	l.push_back(0);
	v.push_back(-1);
	_c.push_back(v.back());
	r = MSK_putbarcblocktriplet(task, (int)v.size(), j.data(), k.data(), l.data(), v.data());
	if (r != MSK_RES_OK)throw std::invalid_argument("MSK_putbarablocktriplet ");

	i.clear();
	j.clear();
	k.clear();
	l.clear();
	v.clear();
	// no constraint for alpha==0
	for (int alphaId(1); alphaId < _id2monomial.size(); ++alphaId) {
		ComplexMonomialPtr const alpha(_id2monomial[alphaId]);
		for (auto const & ij : _B_alpha[alpha]) {
			i.push_back(alphaId - 1);
			j.push_back(id_barvar);
			k.push_back(ij.first.first);
			l.push_back(ij.first.second);
			v.push_back(1);
		}
	}
	r = MSK_putbarablocktriplet(task, (int)v.size(), i.data(), j.data(), k.data(), l.data(), v.data());
	if (r != MSK_RES_OK)throw std::invalid_argument("MSK_putbarablMSK_putbarablocktripletocktriplet ");
	//d::invalid_argument("MSK_putbarablocktriplet ");
	i.clear();
	j.clear();
	k.clear();
	l.clear();
	v.clear();

	i.assign(numcon(), 0);
	v.assign(numcon(), 0);
	for (int n(0); n < i.size(); ++n) {
		i[n] = n;
	}
	std::vector<MSKboundkeye> fx(numcon(), MSK_BK_FX);
	_rhs.assign(numcon(), 0);
	for (auto const & term : _pop->minimize().terms()) {
		//std::cout << "adding " << term.first << std::endl;
		ComplexMonomialPtr2Int::const_iterator ite(_monomial2id.find(term.first));
		if (ite->second != 0) {
			if (ite != _monomial2id.end()) {
				v[ite->second - 1] = term.second.real();
				_rhs[ite->second - 1] = term.second.real();
			}
			else {
				throw std::invalid_argument("unfound monomial ");
			}
		}
	}
	r = MSK_putconboundlist(task, numcon(), i.data(), fx.data(), v.data(), v.data());
	if (r != MSK_RES_OK)throw std::invalid_argument("MSK_putconboundlist ");
}

void SosDualProblem::add_ctr(MSKenv_t & env, MSKtask_t task, int id_ctr) {
	Constraint const & ctr(_pop->ctr(id_ctr));
	double const lb(ctr.lb().real());
	double const ub(ctr.ub().real());
	MSKrescodee  r = MSK_RES_OK;
	if (ub != lb && ub != posInfinity() && lb != negInfinity()) {
		throw std::invalid_argument("RANGE CONSTRAINT DETECTED");
	}

	if (lb == ub) {
		// equality constraint results in free variables
	}
	else {
		// inequality constraint results in sdp variables (if <= constraint is multiplied by -1)
		int id_barvar;
		r = MSK_getnumbarvar(task, &id_barvar);
		if (r != MSK_RES_OK)throw std::invalid_argument("MSK_getnumbarvar ");
		// alpha
		IntVector i;
		// barvar
		IntVector j;
		// Ci_alpha
		IntVector k;
		IntVector l;
		NumberVector v;
		double const factor(ub != posInfinity() ? -1 : 1);
		double const bound(ub != posInfinity() ? ub : lb);
		ComplexPolynomial f_reduced = factor*ctr.f() - bound*factor;


		int const degree(ctr.f().degree());
		int const highest(_monomial2id.rbegin()->first->degree());

		//std::cout << "degree  " << degree << std::endl;
		//std::cout << "highest " << highest << std::endl;
		//std::cout << "factor  " << factor << std::endl;
		//std::cout << "f_reduced = " << f_reduced << std::endl;
		int n(-1);
		for (auto const & term : f_reduced.terms()) {
			int alpha_i = 0;
			for (ComplexMonomialPtr2Int::const_iterator it(_monomial2id.begin()); it != _monomial2id.end() && 2 * it->first->degree() + degree <= highest; ++it, ++alpha_i) {
				int alpha_j = alpha_i;
				for (ComplexMonomialPtr2Int::const_iterator jt(it); jt != _monomial2id.end() && 2 * jt->first->degree() + degree <= highest; ++jt, ++alpha_j) {
					n = std::max(alpha_i, n);
					n = std::max(alpha_j, n);
					ComplexMonomialPtr alpha = it->first + jt->first;
					ComplexMonomialPtr beta  = alpha + term.first;
					//std::cout << "----------------------------" << std::endl;
					//std::cout << it->first << " | " << jt->first <<" : "<<term.second.real() << " "<<beta<< std::endl;
					//std::cout << "alpha_i = " << alpha_i << std::endl;
					//std::cout << "alpha_j = " << alpha_j << std::endl;
					//std::cout << term.second.real() << " " << beta << std::endl;
					if (beta->degree() > 0 && !isZero(term.second.real())) {
						ComplexMonomialPtr2Int::const_iterator alphaIt(_monomial2id.find(alpha));
						ComplexMonomialPtr2Int::const_iterator betaIt(_monomial2id.find(beta));
						if (alphaIt != _monomial2id.end() && betaIt != _monomial2id.end()) {
							C_i_alpha[id_ctr][beta][{alpha_j, alpha_i}] = term.second.real();
						}
						else {
							std::cout << "Unfound  " << alpha << std::endl;
							std::cout << "term     " << term.first << std::endl;
							std::cout << "it       " << it->first << ", " << it->first->degree() << std::endl;
							std::cout << "jt       " << jt->first << ", " << jt->first->degree() << std::endl;
							std::cout << "degree   " << degree << std::endl;
							std::cout << "highest  " << highest << std::endl;
							throw std::invalid_argument("ERROR");
						}
					}
					//_B_alpha[*it->first + *jt->first].insert({ jt->second, it->second });
				}
			}
		}
		_lenbarvar.push_back(n + 1);
		r = MSK_appendbarvars(task, 1, &_lenbarvar.back());
		if (r != MSK_RES_OK)throw std::invalid_argument("MSK_appendbarvars ");
		i.clear();
		j.clear();
		k.clear();
		l.clear();
		v.clear();
		for (auto const & term : C_i_alpha[id_ctr]) {
			ComplexMonomialPtr2Int::const_iterator alphaIt(_monomial2id.find(term.first));
			for (auto const & ij : term.second) {
				i.push_back(alphaIt->second - 1);
				j.push_back(id_barvar);
				k.push_back(ij.first.first);
				l.push_back(ij.first.second);
				v.push_back(ij.second);
			}
		}
		r = MSK_putbarablocktriplet(task, (int)v.size(), i.data(), j.data(), k.data(), l.data(), v.data());
		if (r != MSK_RES_OK)throw std::invalid_argument("MSK_putbarablocktriplet ");
		i.clear();
		j.clear();
		k.clear();
		l.clear();
		v.clear();
		j.push_back(id_barvar);
		k.push_back(0);
		l.push_back(0);
		v.push_back(-f_reduced.constant().real());
		_c.push_back(v.back());
		r = MSK_putbarcblocktriplet(task, (int)v.size(), j.data(), k.data(), l.data(), v.data());
		//std::cout << "C_i_alpha : " << std::endl;
		//for (auto const & term : C_i_alpha[id_ctr]) {
		//	std::cout << std::setw(8) << term.first;
		//	std::cout << " : ";
		//	for (auto const & ij : term.second) {
		//		std::cout << ij.second << "(" << ij.first.first << ", " << ij.first.second << ") ";
		//	}
		//	std::cout << std::endl;
		//}


		//std::exit(0);
	}

}

std::ostream & SosDualProblem::print(std::ostream & stream)const {
	std::cout << "Printing SosDualProblem" << std::endl;
	for (int i(1); i < _id2monomial.size(); ++i) {
		ComplexMonomialPtr alpha(_id2monomial[i]);
		std::cout << std::setw(8) << alpha << std::endl;
		std::cout << std::setw(8) << "RHS";
		std::cout << std::setw(8) << _rhs[i - 1] << std::endl;
		{
			MonomialDecomposition::const_iterator it(_B_alpha.find(alpha));
			if (it != _B_alpha.end()) {
				for (auto const & ij : it->second) {
					std::cout << std::setw(8) << "X";
					std::cout << std::setw(8) << ij.first.first;
					std::cout << std::setw(8) << ij.first.second;
					std::cout << std::setw(8) << ij.second;
					std::cout << std::endl;
				}
			}
		}
		for (int id_ctr(0); id_ctr < C_i_alpha.size(); ++id_ctr) {
			MonomialDecomposition const & ctr(C_i_alpha[id_ctr]);
			MonomialDecomposition::const_iterator it(ctr.find(alpha));
			if (it != ctr.end()) {
				for (auto const & ij : it->second) {
					std::stringstream buffer;
					buffer << "Z_" << id_ctr;
					std::cout << std::setw(8) << buffer.str();
					std::cout << std::setw(8) << ij.first.first;
					std::cout << std::setw(8) << ij.first.second;
					std::cout << std::setw(8) << ij.second;
					std::cout << std::endl;
				}
			}
		}
	}
	for (int i(0); i < numbarvar(); ++i) {
		std::cout << std::setw(8) << "MAX";
		std::stringstream buffer;
		if (i > 0)
			buffer << "Z_" << i - 1;
		else buffer << "X";
		std::cout << std::setw(8) << buffer.str();
		std::cout << std::setw(8) << _c[i];
		std::cout << std::endl;
	}
	return stream;
}


void SosDualProblem::build_obj() {
	IndexedPool2Square & X = _conic.new_block("X", lenbarvar(), true);
	_conic.add_sdp(0, X.id(), 0, 0, -1);
	_c.push_back(-1);
	// no constraint for alpha==0
	for (int alphaId(1); alphaId < _id2monomial.size(); ++alphaId) {
		ComplexMonomialPtr const alpha(_id2monomial[alphaId]);
		for (auto const & ij : _B_alpha[alpha]) {
			_conic.add_sdp(alphaId, X.id(), ij.first.first, ij.first.second, 1);
		}
	}
	for (int n(0); n < numcon(); ++n)
		_conic.new_ctr(0);
	for (auto const & term : _pop->minimize().terms()) {
		//std::cout << "adding " << term.first << std::endl;
		ComplexMonomialPtr2Int::const_iterator ite(_monomial2id.find(term.first));
		if (ite->second != 0) {
			if (ite != _monomial2id.end()) {
				_conic.rhs(ite->second) = term.second.real();
			}
			else {
				throw std::invalid_argument("unfound monomial ");
			}
		}
	}
}
void SosDualProblem::build_ctr(int id_ctr) {
	Constraint const & ctr(_pop->ctr(id_ctr));
	double const lb(ctr.lb().real());
	double const ub(ctr.ub().real());
	if (ub != lb && ub != posInfinity() && lb != negInfinity()) {
		throw std::invalid_argument("RANGE CONSTRAINT DETECTED");
	}

	if (lb == ub) {
		// equality constraint results in free variables
	}
	else {
		// inequality constraint results in sdp variables (if <= constraint is multiplied by -1)
		double const factor(ub != posInfinity() ? -1 : 1);
		double const bound(ub != posInfinity() ? ub : lb);
		ComplexPolynomial f_reduced = factor*ctr.f() - bound*factor;

		int const degree(ctr.f().degree());
		int const highest(_monomial2id.rbegin()->first->degree());

		//std::cout << "degree  " << degree << std::endl;
		//std::cout << "highest " << highest << std::endl;
		//std::cout << "factor  " << factor << std::endl;
		//std::cout << "f_reduced = " << f_reduced << std::endl;
		int n(-1);
		for (auto const & term : f_reduced.terms()) {
			int alpha_i = 0;
			for (ComplexMonomialPtr2Int::const_iterator it(_monomial2id.begin()); it != _monomial2id.end() && 2 * it->first->degree() + degree <= highest; ++it, ++alpha_i) {
				int alpha_j = alpha_i;
				for (ComplexMonomialPtr2Int::const_iterator jt(it); jt != _monomial2id.end() && 2 * jt->first->degree() + degree <= highest; ++jt, ++alpha_j) {
					n = std::max(alpha_i, n);
					n = std::max(alpha_j, n);
					ComplexMonomialPtr alpha = it->first + jt->first;
					ComplexMonomialPtr beta = alpha + term.first;
					//std::cout << "----------------------------" << std::endl;
					//std::cout << it->first << " | " << jt->first <<" : "<<term.second.real() << " "<<beta<< std::endl;
					//std::cout << "alpha_i = " << alpha_i << std::endl;
					//std::cout << "alpha_j = " << alpha_j << std::endl;
					//std::cout << term.second.real() << " " << beta << std::endl;
					if (beta->degree() > 0 && !isZero(term.second.real())) {
						ComplexMonomialPtr2Int::const_iterator alphaIt(_monomial2id.find(alpha));
						ComplexMonomialPtr2Int::const_iterator betaIt(_monomial2id.find(beta));
						if (alphaIt != _monomial2id.end() && betaIt != _monomial2id.end()) {
							C_i_alpha[id_ctr][beta][{alpha_j, alpha_i}] = term.second.real();
						}
						else {
							std::cout << "Unfound  " << alpha << std::endl;
							std::cout << "term     " << term.first << std::endl;
							std::cout << "it       " << it->first << ", " << it->first->degree() << std::endl;
							std::cout << "jt       " << jt->first << ", " << jt->first->degree() << std::endl;
							std::cout << "degree   " << degree << std::endl;
							std::cout << "highest  " << highest << std::endl;
							throw std::invalid_argument("ERROR");
						}
					}
				}
			}
		}
		IndexedPool2Square & Z = _conic.new_block(Str("Z_", id_ctr), n + 1, true);
		for (auto const & term : C_i_alpha[id_ctr]) {
			ComplexMonomialPtr2Int::const_iterator alphaIt(_monomial2id.find(term.first));
			for (auto const & ij : term.second) {
				_conic.add_sdp(alphaIt->second, Z.id(), ij.first.first, ij.first.second, ij.second);
			}
		}
		_conic.add_sdp(0, Z.id(), 0, 0, -f_reduced.constant().real());
	}

}