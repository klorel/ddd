
#include "common.h"
#include "Problem.h"
#include "OpfElement.h"

#include "test_cholesky.h"
#include "RegisteredInstance.h"

#include "MomentGenerator.h"
#include "SdpProblem.h"
#include "RealPolynomial.h"

int test_pooling(int, char**);
int test_sos(int, char**);
int test_sdp(int, char**);

void get_sdp_1(SdpProblem & sdp) {
	sdp.clear();
	int const b1 = sdp.newBlock(2);
	int const b2 = sdp.newBlock(3);
	int const b3 = sdp.newBlock(-2);

	int y1 = sdp.newy(1);
	int y2 = sdp.newy(2);

	sdp.add(0, b1, 1, 1, 2);
	sdp.add(0, b1, 1, 2, 1);
	sdp.add(0, b1, 2, 2, 2);

	sdp.add(0, b2, 1, 1, 3);
	sdp.add(0, b2, 1, 3, 1);
	sdp.add(0, b2, 2, 2, 2);
	sdp.add(0, b2, 3, 3, 3);


	sdp.add(y1, b1, 1, 1, 3);
	sdp.add(y1, b1, 1, 2, 1);
	sdp.add(y1, b1, 2, 2, 3);

	sdp.add(y1, b3, 1, 1, 1);

	sdp.add(y2, b2, 1, 1, 3);
	sdp.add(y2, b2, 1, 3, 1);
	sdp.add(y2, b2, 2, 2, 4);
	sdp.add(y2, b2, 3, 3, 5);

	sdp.add(y2, b3, 2, 2, 1);

	int t = sdp.newBlock(-2);
	sdp.add(1, t, 1, 1, -1);
	sdp.add(1, t, 2, 2, +1);
	sdp.add(2, t, 1, 1, -1);
	sdp.add(2, t, 2, 2, +1);

}
// attention Qk et 0.5 Lk
void get_sdp_2(SdpProblem & sdp) {
	sdp.clear();
	sdp.newBlock(2);
	sdp.add(0, 1, 1, 2, +1);
	sdp.add(0, 1, 2, 2, -1);
	sdp.newBlock(-1);
	sdp.add(0, 2, 1, 1, -1);
	sdp.newy(1);
	sdp.add(1, 1, 1, 1, 1);
}
void get_sdp_3(SdpProblem & sdp) {
	sdp.clear();
	sdp.newBlock(3);
	// max - X11 
	sdp.add(0, 1, 1, 1, -1);
	// 0 = X12 = y10
	sdp.newy(0);
	sdp.add(1, 1, 1, 3, 1);
	// -2 = X13 = y01
	sdp.newy(-2);
	sdp.add(2, 1, 1, 2, 1);
	// +0 = X23 = y11
	sdp.newy(0);
	sdp.add(3, 1, 2, 3, 1);
	// +1 = X22 = y20
	sdp.newy(+1);
	sdp.add(4, 1, 3, 3, 1);
	// +1 = X33 = y02
	sdp.newy(+1);
	sdp.add(5, 1, 2, 2, 1);

}
int main(int argc, char** argv) {
	//SdpProblem sdp;
	//get_sdp_2(sdp);
	//sdp.print("toto.dat");

	test_pooling(argc, argv);
	return 0;
}

int test_sdp(int argc, char**) {

	return 0;
}
int test_sos(int argc, char** argv) {

	int n(2);
	int d(2);
	if (argc >= 2) {
		std::stringstream buffer(argv[1]);
		buffer >> n;
	}
	if (argc >= 3) {
		std::stringstream buffer(argv[2]);
		buffer >> d;
	}
	MomentGenerator generator(n, d);
	generator.build();
	return 0;
}
int test_pooling(int argc, char**argv) {
	Problem pop;
	pop.newvarpool("x", 1);
	pop.newvarpool("y", 1);

	RealPolynomial x;
	x.set(0);
	RealPolynomial y;
	y.set(1);

	RealPolynomial m;

	m = x;
	m = m + y;
	m = m + x*y;
	m = m + y*y;
	m = m + y*x;
	m.print(std::cout, pop);

	return 0;
	AvailableInstances id(AvailableInstances::SIZE);
	std::cout << "argc = " << argc << std::endl;
	if (argc >= 2) {
		std::stringstream buffer(argv[1]);
		size_t i;
		buffer >> i;
		id = static_cast<AvailableInstances>(i);
	}
	int useP(false);
	if (argc >= 3) {
		std::stringstream buffer(argv[2]);
		buffer >> useP;
	}
	if (id == AvailableInstances::SIZE) {
		RegisteredInstance::PrintAvailable();
		return 0;
	}
	RegisteredInstance instance(id);
	Problem problem;
	if (useP)
		instance.pFormulation(problem);
	else
		instance.pqFormulation(problem);
	problem.amplExport("problem");
	return 0;
	//std::cout << "pq formulation : " << std::endl << problem << std::endl;
	//for (size_t n(0); n < problem.nvars(); ++n){
	//	if (!sparsityPattern[n].empty()){
	//		std::cout << problem.name(n) << " : ";
	//		for (auto const & m : sparsityPattern[n])
	//			std::cout << problem.name(m) << " ";
	//		std::cout << std::endl;
	//	}
	//}
	std::cout << "problem is "<<std::endl<<problem << std::endl;
	problem.removeInequality();
	std::cout << "removed inequality" << std::endl<< problem << std::endl;

	SdpProblem sdp;
	sdp.sdprelaxation(problem);
	sdp.print("sdprelaxation.dat");
	sdp.print(std::cout, "X");
	SdpProblem ssdp;
	sdp.sparsesdp(ssdp);
	ssdp.print("ssdp.dat");
	//ssdp.print(std::cout, "Y");
	return 0;
}

//std::cout << "total_size is " << total_size << std::endl;
//if (total_size < 50) {
//	for (auto const & clique : cliqueDecomposition) {
//		for (auto const & i : *clique) {
//			if(i>0)
//			std::cout << problem.name(i-1) << " ";
//		}
//		std::cout << std::endl;
//	}
//}
//std::cout << "matrix.size() is " << matrix.cols() << std::endl;
//std::cout << "max_size      is " << max_size << std::endl;
//std::cout << "nb clique     is " << cliqueDecomposition.size() << std::endl;
//size_t max_support(0);
//SparsityPattern support;
//problem.addSupport(support);
//for (int i(0); i < problem.nctrs() + 1; ++i) {
//	if (support[i].size() > max_support) {
//		max_support = support[i].size();
//		if (i < problem.nctrs())
//			std::cout << "max support becomes " << max_support << ", ctr[" << i << "]:" << problem.ctrname(i) << std::endl;
//		else
//			std::cout << "max support becomes " << max_support << ", obj " << std::endl;
//	}
//}
//std::cout << "max_support   is " << max_support << std::endl;