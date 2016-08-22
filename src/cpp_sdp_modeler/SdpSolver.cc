#include "SdpSolver.h"

#include <xprs.h>

#include "mosek.h"    /* Include the MOSEK definition file.  */

#include <SparseCore>
#include <SparseCholesky>
#include <Eigenvalues>


typedef Eigen::SparseMatrix<double, 0, int> SparseMatrix;
typedef Eigen::Triplet<double> Triplet;
typedef std::vector<Triplet> Triplets;


SdpSolver::SdpSolver(SdpProblem const & input) :_input(input) {

}

SdpSolver::~SdpSolver() {

}


static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
	std::cout << str;;
}


void SdpSolver::launch_mosek() {
	int const NUMCON = (int)_input._b.size();
	int const NUMBARVAR = (int)_input._blocks.size();
	std::vector<int> DIMBARVAR(NUMBARVAR);
	std::vector<int> LENBARVAR(NUMBARVAR);
	for (size_t i(0); i < _input._blocks.size(); ++i) {
		auto const & block(_input._blocks[i]);
		DIMBARVAR[i] = std::abs(block._size);
		LENBARVAR[i] = block._size > 0 ? DIMBARVAR[i] * (DIMBARVAR[i] + 1) / 2 : DIMBARVAR[i];
	}

	//
	MSKrescodee  r;

	MSKenv_t     env = NULL;
	MSKtask_t    task = NULL;
	/* Create the mosek environment. */
	r = MSK_makeenv(&env, NULL);
	if (r == MSK_RES_OK)
	{
		/* Create the optimization task. */
		r = MSK_maketask(env, NUMCON, 0, &task);
		if (r == MSK_RES_OK)
		{
			MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);

			/* Append 'NUMCON' empty constraints.
			The constraints will initially have no bounds. */
			if (r == MSK_RES_OK)
				r = MSK_appendcons(task, NUMCON);
			/* Append 'NUMBARVAR' semidefinite variables. */
			if (r == MSK_RES_OK) {
				r = MSK_appendbarvars(task, NUMBARVAR, DIMBARVAR.data());
			}
			// sorted by ctr then by block and then UPPER triangle 1-based
			{
				IntVector bara_i;
				IntVector bara_j;
				IntVector bara_k;
				IntVector bara_l;
				NumberVector bara_v;

				IntVector barc_j;
				IntVector barc_k;
				IntVector barc_l;
				NumberVector barc_v;
				for (auto const & kvp : _input._matrix) {
					int const ctr(kvp.first[0]);
					int const block(kvp.first[1]);
					int const xi(kvp.first[2]);
					int const xj(kvp.first[3]);
					if (ctr == 0) {
						barc_j.push_back(block - 1);
						barc_k.push_back(xj - 1);
						barc_l.push_back(xi - 1);
						barc_v.push_back(kvp.second);
					}
					else {
						bara_i.push_back(ctr - 1);
						bara_j.push_back(block - 1);
						bara_k.push_back(xj - 1);
						bara_l.push_back(xi - 1);
						bara_v.push_back(kvp.second);
					}
				}
				std::cout << "MSK_putbarcblocktriplet " << std::endl;
				if (r == MSK_RES_OK)
					r = MSK_putbarcblocktriplet(task, (int)barc_v.size(), barc_j.data(), barc_k.data(), barc_l.data(), barc_v.data());
				std::cout << "MSK_putbarablocktriplet " << std::endl;
				if (r == MSK_RES_OK)
					r = MSK_putbarablocktriplet(task, (int)bara_v.size(), bara_i.data(), bara_j.data(), bara_k.data(), bara_l.data(), bara_v.data());
			}
		}
		if (r == MSK_RES_OK)
			r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);
		{
			//MSK_putconbound
			std::vector<MSKboundkey_enum> fx(NUMCON, MSK_BK_FX);
			IntVector conidx(NUMCON);
			NumberVector bnd(NUMCON);
			for (int i(0); i < NUMCON; ++i) {
				conidx[i] = i;
				bnd[i] = _input._b[i];
			}
			std::cout << "MSK_putconboundlist " << std::endl;
			r = MSK_putconboundlist(task, NUMCON, conidx.data(), fx.data(), bnd.data(), bnd.data());
		}

		if (r == MSK_RES_OK)
		{
			MSKrescodee trmcode;

			/* Run optimizer */
			r = MSK_optimizetrm(task, &trmcode);

			/* Print a summary containing information
			about the solution for debugging purposes*/
			MSK_solutionsummary(task, MSK_STREAM_MSG);

			NumberVector barx;

			if (r == MSK_RES_OK)
			{
				MSKsolstae solsta;

				MSK_getsolsta(task, MSK_SOL_ITR, &solsta);

				switch (solsta)
				{
				case MSK_SOL_STA_OPTIMAL:
				case MSK_SOL_STA_NEAR_OPTIMAL:
					for (int i(0); i < NUMBARVAR; ++i) {
						barx.assign(LENBARVAR[i], -1);
						MSK_getbarxj(task,
							MSK_SOL_ITR,    /* Request the interior solution. */
							i,
							barx.data());
						//for (int j(0); j < LENBARVAR[i]; ++j) {
						//	std::cout << std::setw(5) << i;
						//	std::cout << std::setw(5) << j;
						//	std::cout << std::setw(15) << barx[j] << std::endl;;
						//}
					}


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
		/* Delete the task and the associated data. */
		MSK_deletetask(&task);
	}
	/* Delete the environment and the associated data. */
	MSK_deleteenv(&env);

}


void XPRS_CC optimizermsg(XPRSprob prob, void* data, const char *sMsg, int nLen, int nMsgLvl);
void errormsg(XPRSprob const & prob, const char *sSubName, int nLineNo, int nErrorCode);

bool is_sdp_set(IntVector const & sdp_set, StrVector const & col_name, NumberVector const & x, IntVector & cut_start, IntVector & cut_index, NumberVector & cut_value, size_t max_cut);

void SdpSolver::launch_xpress() {
	int nReturn;
	char banner[256];
	//
	XPRSinit(NULL);
	XPRSgetbanner(banner);
	std::cout << banner << std::endl;
	XPRSprob prob;
	nReturn = XPRScreateprob(&prob);
	if (nReturn != 0 && nReturn != 32) errormsg(prob, "XPRScreateprob", __LINE__, nReturn);
	// creating empty problem
	nReturn = XPRSloadlp(prob, "sdp_cutting", 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
	if (nReturn != 0)errormsg(prob, "XPRSloadlp", __LINE__, nReturn);
	int ncols(_input.nz());
	int nrows((int)_input._b.size());
	// adding columns

	std::cout << "ncols is " << ncols << std::endl;
	// X,Z,y
	NumberVector col_lb(ncols * 2 + nrows, -1e6);
	NumberVector col_ub(col_lb.size(), +1e6);
	NumberVector col_obj(col_lb.size(), 0);
	std::cout << "nrows is " << nrows << std::endl;
	// AX=A, Z=yA-C
	NumberVector row_rhs(nrows + ncols);
	std::copy(_input._b.begin(), _input._b.end(), row_rhs.begin());
	std::vector<char> row_sense(row_rhs.size(), 'E');

	IntVector row_start;
	IntVector col_index;
	NumberVector mat_value;

	// PRIMAL PART, AX=A
	for (auto const & kvp : _input._matrix) {
		int const ctr(kvp.first[0]);
		int const i(kvp.first[1]);
		int const j(kvp.first[2]);
		int const k(kvp.first[3]);
		//
		int const id(_input.id(i - 1, j - 1, k - 1));
		if (j == k) {
			col_lb[id] = 0;
			col_lb[ncols + id] = 0;
		}
		if (ctr == 0) {
			col_obj[id] = kvp.second*0;
		}
		else {
			if (row_start.size() != ctr) {
				// new row
				row_start.push_back((int)mat_value.size());
			}
			col_index.push_back(id);
			mat_value.push_back(kvp.second);
		}
	}
	Matrix dual;
	_input.dual(dual);
	int last_zid(-1);
	// DUAL PART, Z - yA = -C
	for (auto const & kvp : dual) {
		int const b(kvp.first[0]);
		int const bi(kvp.first[1]);
		int const bj(kvp.first[2]);
		int const ctr(kvp.first[3]);
		//
		int const id(_input.id(b - 1, bi - 1, bj - 1));
		int const zid(ncols + id);
		// new row
		if (last_zid != zid) {
			row_start.push_back((int)mat_value.size());
			col_index.push_back(zid);
			mat_value.push_back(1);
			last_zid = zid;
		}
		if (ctr == 0) {
			row_rhs[row_start.size() - 1] = -kvp.second;
		}
		else {
			//
			int const yid(2 * ncols + ctr - 1);
			col_index.push_back(yid);
			mat_value.push_back(-kvp.second);
		}
	}
	// weak duality
	row_start.push_back((int)mat_value.size());
	row_rhs.push_back(0);
	row_sense.push_back('G');
	for (int i(0); i < ncols; ++i) {
		col_index.push_back(i);
		mat_value.push_back(-1);
	}
	for (int i(0); i < nrows; ++i) {
		col_index.push_back(2 * ncols + i);
		mat_value.push_back(row_rhs[i]);
		col_obj[2 * ncols + i] = row_rhs[i];
	}
	//
	nReturn = XPRSaddcols(prob, (int)col_lb.size(), 0, col_obj.data(), NULL, NULL, NULL, col_lb.data(), col_ub.data());
	if (nReturn != 0)errormsg(prob, "XPRSaddcols", __LINE__, nReturn);
	//
	row_start.push_back((int)mat_value.size());
	nReturn = XPRSaddrows(prob, (int)row_rhs.size(), (int)mat_value.size(), row_sense.data(), row_rhs.data(), NULL, row_start.data(), col_index.data(), mat_value.data());
	if (nReturn != 0)errormsg(prob, "XPRSaddrows", __LINE__, nReturn);


	// sdp sets
	std::vector<IntVector> sdp_cones;
	for (int i(0); i < _input.nblock(); ++i) {
		SdpProblem::Block const & block(_input._blocks[i]);
		if (block._size > 0) {
			IntVector primal_set;
			IntVector dual_set;
			for (int j(0); j < block._size; ++j)
				for (int k(j); k < block._size; ++k) {
					int const id(_input.id(i, j, k));
					primal_set.push_back(id);
					dual_set.push_back(id + ncols);
				}
			//sdp_cones.push_back(primal_set);
			sdp_cones.push_back(dual_set);
		}
	}

	StrVector col_name(col_lb.size());
	for (int i(0); i < _input.nblock(); ++i) {
		SdpProblem::Block const & block(_input._blocks[i]);
		if (block._size > 0) {
			for (int j(0); j < block._size; ++j)
				for (int k(j); k < block._size; ++k) {
					int const id(_input.id(i, j, k));
					//std::cout << "id = " << id << std::endl;
					col_name[id] = Str("X_", i, "_", j, "_", k);
					col_name[ncols + id] = Str("Z_", i, "_", j, "_", k);
				}
		}
		else {
			for (int j(0); j < -block._size; ++j) {
				int const id(_input.id(i, j, j));
				col_name[id] = Str("X_", i, "_", j, "_", j);
				col_name[ncols + id] = Str("Z_", i, "_", j, "_", j);
			}
		}
	}
	for (int i(0); i < nrows; ++i) {
		col_name[2 * ncols + i] = Str("y_", i);
	}
	for (int i(0); i < col_lb.size(); ++i) {
		//std::cout << "col_name["<<i<<"] = "<<col_name[i] << std::endl;
		nReturn = XPRSaddnames(prob, 2, col_name[i].c_str(), i, i);
		if (nReturn != 0)
			errormsg(prob, "XPRSaddnames", __LINE__, nReturn);
	}


	nReturn = XPRSsetintcontrol(prob, XPRS_DEFAULTALG, XPRS_ALG_DUAL);
	if (nReturn != 0)errormsg(prob, "XPRSsetintcontrol", __LINE__, nReturn);

	nReturn = XPRSsetintcontrol(prob, XPRS_THREADS, 1);
	if (nReturn != 0)errormsg(prob, "XPRSsetintcontrol", __LINE__, nReturn);

	nReturn = XPRSsetcbmessage(prob, optimizermsg, NULL);
	if (nReturn != 0)errormsg(prob, "XPRSsetcbmessage", __LINE__, nReturn);

	nReturn = XPRSchgobjsense(prob, XPRS_OBJ_MINIMIZE);
	if (nReturn != 0)errormsg(prob, "XPRSchgobjsense", __LINE__, nReturn);

	//nReturn = XPRSsetdblcontrol(prob, XPRS_FEASTOL, 1e-8);
	//if (nReturn != 0)errormsg(prob, "XPRSsetdblcontrol", __LINE__, nReturn);
	//nReturn = XPRSsetdblcontrol(prob, XPRS_FEASTOLTARGET, 1e-8);
	//if (nReturn != 0)errormsg(prob, "XPRSsetdblcontrol", __LINE__, nReturn);
	//nReturn = XPRSsetdblcontrol(prob, XPRS_OPTIMALITYTOL, 0);
	//if (nReturn != 0)errormsg(prob, "XPRSsetdblcontrol", __LINE__, nReturn);
	//nReturn = XPRSsetdblcontrol(prob, XPRS_OPTIMALITYTOLTARGET, 0);
	//if (nReturn != 0)errormsg(prob, "XPRSsetdblcontrol", __LINE__, nReturn);

	nReturn = XPRSsetintcontrol(prob, XPRS_SCALING, 0);
	if (nReturn != 0)errormsg(prob, "XPRSsetintcontrol", __LINE__, nReturn);
	//nReturn = XPRSsetintcontrol(prob, XPRS_PRESOLVE, 0);
	//if (nReturn != 0)errormsg(prob, "XPRSsetintcontrol", __LINE__, nReturn);
	//XPRSwriteprob(prob, "sdp_cutting.lp", "lp");

	XPRSsetintcontrol(prob, XPRS_LPLOG, XPRS_OUTPUTLOG_NO_OUTPUT);

	size_t ite(0);
	bool stop(false);
	NumberVector x(col_lb.size());
	while (!stop) {
		++ite;
		XPRSlpoptimize(prob, "");
		XPRSgetlpsol(prob, x.data(), NULL, NULL, NULL);

		IntVector cut_start;
		IntVector cut_index;
		NumberVector cut_value;
		size_t max_cut(50);
		for (auto const & cone: sdp_cones) {
			is_sdp_set(cone, col_name, x, cut_start, cut_index, cut_value, max_cut);
			if (max_cut != 0 && cut_start.size() >= max_cut)
				break;
		}
		NumberVector cut_rhs(cut_start.size(), 0);
		std::vector<char> cut_sense(cut_start.size(), 'G');
		if (!cut_rhs.empty()) {
			//std::cout << "Number of cuts : " << cut_rhs.size()<<", elements "<<cut_value.size()<<", density "<< cut_value.size()*1.0/cut_rhs.size() << std::endl;
			cut_start.push_back((int)cut_value.size());
			nReturn = XPRSaddrows(prob, (int)cut_rhs.size(), (int)cut_value.size(), cut_sense.data(), cut_rhs.data(), NULL, cut_start.data(), cut_index.data(), cut_value.data());
		}
		//for (int i(0); i < _input.nblock(); ++i) {
		//	SdpProblem::Block const & block(_input._blocks[i]);
		//	if (block._size > 0) {
		//		std::cout << "block._size is " << block._size << std::endl;
		//		int id(block._begin);
		//		Triplets triplets;
		//		for (int j(0); j < block._size; ++j) {
		//			for (int k(j); k < block._size; ++k, ++id) {
		//				if (!isZero(x[id])) {
		//					std::cout << col_name[id] << " : " << x[id] << std::endl;
		//					triplets.push_back({ j,k,x[id] });
		//					if (k != j)
		//						triplets.push_back({ k,j,x[id] });
		//				}
		//			}
		//		}
		//		Eigen::SparseMatrix<double> m(block._size, block._size);
		//		m.setFromTriplets(triplets.begin(), triplets.end());
		//		Eigen::MatrixXd a(m);
		//		Eigen::EigenSolver<Eigen::MatrixXd> eigenof(a, true);
		//		Eigen::EigenSolver<Eigen::MatrixXd>::EigenvalueType eigenvalues = eigenof.eigenvalues();
		//		Eigen::EigenSolver<Eigen::MatrixXd>::EigenvectorsType eigenvectors = eigenof.eigenvectors();
		//		std::cout << "A is " << std::endl << a << std::endl;
		//		for (int lambdaIdx(0); lambdaIdx < eigenvectors.cols(); ++lambdaIdx) {
		//			double const lambda(eigenvalues(lambdaIdx).real());
		//			if (lambda < -1e-10) {
		//				std::cout << "lambda[" << std::setw(6) << lambdaIdx << "] = ";
		//				std::cout << std::setw(15) << lambda;
		//				std::cout << std::endl;
		//				for (int vectorIdx(0); vectorIdx < eigenvectors.rows(); ++vectorIdx) {
		//					std::cout << "v[" << std::setw(6) << vectorIdx << "] = " << eigenvectors(vectorIdx, lambdaIdx) << std::endl;
		//				}
		//				std::cout << "eigenvectors.col(" << lambdaIdx << ")" << std::endl << eigenvectors.col(lambdaIdx) << std::endl;
		//				std::cout << "eigenvectors" << std::endl << eigenvectors << std::endl;
		//				cut_rhs.push_back(0);
		//				cut_sense.push_back('G');
		//				cut_start.push_back((int)cut_value.size());
		//				id = block._begin;
		//				double verif_value(0);
		//				double const scale_factor(1 + 0 * 1.0 / std::sqrt(-lambda));
		//				for (int j(0); j < block._size; ++j) {
		//					Number const vj(eigenvectors(j, lambdaIdx).real());
		//					for (int k(j); k < block._size; ++k, ++id) {
		//						Number const factor(j == k ? 1 : 2);
		//						Number const vk(eigenvectors(k, lambdaIdx).real());
		//						cut_index.push_back(id);
		//						cut_value.push_back(factor*vj*vk*scale_factor);
		//						verif_value += cut_value.back()*x[id];
		//						std::cout << col_name[id] << " | " << cut_value.back() << std::endl;
		//					}
		//				}
		//				double  eigen_verif = (eigenvectors.col(lambdaIdx).transpose()*a*eigenvectors.col(lambdaIdx)).real()(0, 0)*scale_factor;
		//				if (std::abs(eigen_verif - verif_value) > 1e-6) {
		//					std::cout << "error eigen_verif " << std::setprecision(15) << eigen_verif << std::endl;
		//					std::cout << "error verif_value " << std::setprecision(15) << verif_value << std::endl;
		//					std::exit(0);
		//				}

		//			}
		//		}
		//	}
		//}
		//if (!cut_rhs.empty()) {
		//	std::cout << "Number of cuts : " << cut_rhs.size() << std::endl;
		//	nReturn = XPRSaddrows(prob, (int)cut_rhs.size(), (int)cut_value.size(), cut_sense.data(), cut_rhs.data(), NULL, cut_start.data(), cut_index.data(), cut_value.data());
		//}
		//XPRSwriteprob(prob, Str("prob_", ite, ".lp").c_str(), "lp");
		//if (nReturn != 0)errormsg(prob, "XPRSwriteprob", __LINE__, nReturn);
		//XPRSwriteprob(prob, Str("sdp_cutting.lp").c_str(), "lp");
		//if (nReturn != 0)errormsg(prob, "XPRSwriteprob", __LINE__, nReturn);
		stop = cut_rhs.empty();
		int simplex_ite;
		double obj_value;
		XPRSgetintattrib(prob, XPRS_SIMPLEXITER, &simplex_ite);
		XPRSgetdblattrib(prob, XPRS_LPOBJVAL, &obj_value);
		std::cout << std::setw(6) << ite;
		std::cout << std::setw(6) << cut_rhs.size();
		std::cout << std::setw(8) << simplex_ite;
		std::cout << std::setw(25)<<std::setprecision(15) << obj_value;
		std::cout << std::endl;
	}
	//nReturn = XPRSwriteprob(prob, "sdp_cutting.lp", "lp");
	//if (nReturn != 0)errormsg(prob, "XPRSwriteprob", __LINE__, nReturn);

	XPRSdestroyprob(prob);
	XPRSfree();
}

bool is_sdp_set(IntVector const & sdp_set, StrVector const & col_name, NumberVector const & x, IntVector & cut_start, IntVector & cut_index, NumberVector & cut_value, size_t max_cut) {
	bool active_log(false);
	int const n((int)std::round(std::sqrt(2 * sdp_set.size() + 0.75) - 0.5));
	if (active_log) {
		std::cout << "s is " << sdp_set.size() << std::endl;
		std::cout << "n is " << n << std::endl;
		std::cout << "sdt_set is ";
		for (int i(0); i < n; ++i) {
			std::cout << col_name[sdp_set.front() + i] << ", ";
		}
		std::cout << std::endl;
	}
	Eigen::MatrixXd a(n, n);
	a.setZero();
	int id(sdp_set.front());
	for (int i(0); i < n; ++i) {
		for (int j(i); j < n; ++j, ++id) {
			a(i, j) = x[id];
			if (i != j)
				a(j, i) = x[id];
		}
	}

	Eigen::EigenSolver<Eigen::MatrixXd> eigenof(a, true);
	Eigen::EigenSolver<Eigen::MatrixXd>::EigenvalueType eigenvalues = eigenof.eigenvalues();
	Eigen::EigenSolver<Eigen::MatrixXd>::EigenvectorsType eigenvectors = eigenof.eigenvectors();
	if (active_log) {
		std::cout << "A is " << std::endl << a << std::endl;
	}
	for (int lambdaIdx(0); lambdaIdx < eigenvectors.cols(); ++lambdaIdx) {
		double const lambda(eigenvalues(lambdaIdx).real());
		if (lambda < -1e-10) {
			if (active_log) {
				std::cout << "lambda[" << std::setw(6) << lambdaIdx << "] = ";
				std::cout << std::setw(15) << lambda;
				std::cout << std::endl;
				for (int vectorIdx(0); vectorIdx < eigenvectors.rows(); ++vectorIdx) {
					std::cout << "v[" << std::setw(6) << vectorIdx << "] = " << eigenvectors(vectorIdx, lambdaIdx) << std::endl;
				}
				std::cout << "eigenvectors.col(" << lambdaIdx << ")" << std::endl << eigenvectors.col(lambdaIdx) << std::endl;
				std::cout << "eigenvectors" << std::endl << eigenvectors << std::endl;
			}
			cut_start.push_back((int)cut_value.size());

			double verif_value(0);
			double const scale_factor(1 + 0 * 1.0 / std::sqrt(-lambda));


			int id = sdp_set.front();
			for (int i(0); i < n; ++i) {
				Number const vi(eigenvectors(i, lambdaIdx).real());
				for (int j(i); j < n; ++j, ++id) {
					Number const vj(eigenvectors(j, lambdaIdx).real());
					Number const factor(i == j ? 1 : 2);
					cut_index.push_back(id);
					cut_value.push_back(factor*vi*vj*scale_factor);
					verif_value += cut_value.back()*x[id];
					if (active_log) {
						std::cout << col_name[id] << " | " << cut_value.back() << std::endl;
					}
				}
			}

			double  eigen_verif = (eigenvectors.col(lambdaIdx).transpose()*a*eigenvectors.col(lambdaIdx)).real()(0, 0)*scale_factor;
			if (std::abs(eigen_verif - verif_value) > 1e-6) {
				std::cout << "error eigen_verif " << std::setprecision(15) << eigen_verif << std::endl;
				std::cout << "error verif_value " << std::setprecision(15) << verif_value << std::endl;
				std::exit(0);
			}

		}
	}
	return false;
}

/**********************************************************************************\
* Name:         optimizermsg                                                           *
* Purpose:      Display Optimizer error messages and warnings.                         *
* Arguments:    const char *sMsg       Message string                                  *
*               int nLen               Message length                                  *
*               int nMsgLvl            Message type                                    *
* Return Value: None                                                                   *
\**********************************************************************************/
void XPRS_CC optimizermsg(XPRSprob prob, void* data, const char *sMsg, int nLen, int nMsgLvl)
{
	switch (nMsgLvl) {

		/* Print Optimizer error messages and warnings */
	case 4:       /* error */
	case 3:       /* warning */
	case 2:       /* dialogue */
	case 1:       /* information */
		//std::cout << sMsg << std::endl;
		break;
	default:
		fflush(NULL);
		break;
	}
}


/**************************************************************************************\
* Name:         errormsg                                                               *
* Purpose:      Display error information about failed subroutines.                    *
* Arguments:    const char *sSubName   Subroutine name                                 *
*               int nLineNo            Line number                                     *
*               int nErrCode           Error code                                      *
* Return Value: None                                                                   *
\**************************************************************************************/
void errormsg(XPRSprob const & prob, const char *sSubName, int nLineNo, int nErrCode)
{
	int nErrNo;   /* Optimizer error number */

				  /* Print error message */
	printf("The subroutine %s has failed on line %d\n", sSubName, nLineNo);

	/* Append the error code, if it exists */
	if (nErrCode != -1) printf("with error code %d.\n\n", nErrCode);

	/* Append Optimizer error number, if available */
	if (nErrCode == 32) {
		XPRSgetintattrib(prob, XPRS_ERRORCODE, &nErrNo);
		printf("The Optimizer error number is: %d.\n\n", nErrNo);
	}

	/* Free memory, close files and exit */
	XPRSdestroyprob(prob);
	XPRSfree();
	exit(nErrCode);
}



