
#include "common.h"
#include "RegisteredInstance.h"
#include "Problem.h"

int main(int argc, char**argv) {
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
	std::cout << problem << std::endl;
	return 0;
	//problem.amplExport("problem");
	////std::cout << "pq formulation : " << std::endl << problem << std::endl;
	////for (size_t n(0); n < problem.nvars(); ++n){
	////	if (!sparsityPattern[n].empty()){
	////		std::cout << problem.name(n) << " : ";
	////		for (auto const & m : sparsityPattern[n])
	////			std::cout << problem.name(m) << " ";
	////		std::cout << std::endl;
	////	}
	////}
	//std::cout << "problem is " << std::endl << problem << std::endl;
	//problem.removeInequality();
	//std::cout << "removed inequality" << std::endl << problem << std::endl;

	//SdpProblem sdp;
	//sdp.sdprelaxation(problem);
	//sdp.print("sdprelaxation.dat");
	//sdp.print(std::cout, "X");
	//SdpProblem ssdp;
	//sdp.sparsesdp(ssdp);
	//ssdp.print("ssdp.dat");
	////ssdp.print(std::cout, "Y");
	return 0;
}