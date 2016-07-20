
#include "common.h"
#include "SdpProblem.h"
#include "SdpSolver.h"

int main(int argc, char**argv) {
	SdpProblem sdp1;
	get_sdp_1(sdp1);
	SdpSolver solver(sdp1);
	//solver.launch_mosek();
	//sdp1.print("sdp1.dat");

	solver.launch_xpress();
}