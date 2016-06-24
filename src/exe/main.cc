
#include "common.h"
#include "Problem.h"
#include "OpfElement.h"

#include "test_cholesky.h"
#include "RegisteredInstance.h"

#include "MomentGenerator.h"

int main(int argc, char** argv) {
	//int n(2);
	//int d(2);
	//if (argc >= 2){
	//	std::stringstream buffer(argv[1]);
	//	buffer >> n;
	//}
	//if (argc >= 3){
	//	std::stringstream buffer(argv[2]);
	//	buffer >> d;
	//}
	//MomentGenerator generator(n,d);
	//generator.build();
	//return 0;

	AvailableInstances id(AvailableInstances::SIZE);
	std::cout << "argc = " << argc << std::endl;
	if (argc >= 2){
		std::stringstream buffer(argv[1]);
		size_t i;
		buffer >> i;
		id = static_cast<AvailableInstances>(i);
	}
	int useP(false);
	if (argc >= 3){
		std::stringstream buffer(argv[2]);
		buffer >> useP;
	}
	if (id == AvailableInstances::SIZE){
		RegisteredInstance::PrintAvailable();
		return 0;
	}
	RegisteredInstance instance(id);
	Problem problem;
	if (useP)
		instance.pFormulation(problem);
	else
		instance.pqFormulation(problem);
	//std::cout << "pq formulation : " << std::endl << problem << std::endl;
	SparsityPattern sparsityPattern;
	problem.addSparsityPattern(sparsityPattern);
	//for (size_t n(0); n < problem.nvars(); ++n){
	//	if (!sparsityPattern[n].empty()){
	//		std::cout << problem.name(n) << " : ";
	//		for (auto const & m : sparsityPattern[n])
	//			std::cout << problem.name(m) << " ";
	//		std::cout << std::endl;
	//	}
	//}
	//std::cout << problem << std::endl;

	SparseMatrix matrix;
	build(sparsityPattern, matrix);

	//std::cout << matrix << std::endl;
	IntSetPtrSet cliqueDecomposition;
	build(matrix, cliqueDecomposition);

	size_t total_size(0);
	size_t max_size(0);
	for (auto const & clique : cliqueDecomposition){
		total_size += clique->size();
		max_size = std::max(max_size, clique->size());
	}
	std::cout << "total_size is " << total_size << std::endl;
	if (total_size < 50){
		for (auto const & clique : cliqueDecomposition){
			for (auto const & i : *clique){
				std::cout << problem.name(i) << " ";
			}
			std::cout << std::endl;
		}
	}
	std::cout << "matrix.size() is " << matrix.cols() << std::endl;
	std::cout << "max_size      is " << max_size << std::endl;
	std::cout << "nb clique     is " << cliqueDecomposition.size() << std::endl;
	size_t max_support(0);
	SparsityPattern support;
	problem.addSupport(support);
	for (int i(0); i<problem.nctrs()+1; ++i){
		if (support[i].size()>max_support){
			max_support = support[i].size();
			if (i<problem.nctrs())
				std::cout << "max support becomes " << max_support << ", ctr["<<i<<"]:" << problem.ctrname(i) << std::endl;
			else
				std::cout << "max support becomes " << max_support << ", obj " << std::endl;
		}
	}
	std::cout << "max_support   is " << max_support << std::endl;
	//test_cholesky(argc, argv);

	//Eigen::NaturalOrdering
	//ElementT<branch> b;
	//Problem p;
	//size_t n(10);
	//p.ivariablePool("v", n);

	// 1+x0+x2
	//FunctionReal f(-1 + p.variable("v_real", 7) - p.variable("v_real", 2)*p.variable("v_imag", 2));
	//f.print(std::cout, p);

	//std::vector<FunctionComplex> allV(n);
	//for (size_t i(0); i < n; ++i){
	//	allV[i] = p.ivariable("v", i);
	//}

	//for (auto & v : allV){
	//	v.print(std::cout, p);
	//	std::cout << std::endl;
	//}
	return 0;
}