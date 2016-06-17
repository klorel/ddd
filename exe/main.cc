
#include "common.h"
#include "Problem.h"
#include "OpfElement.h"

int main(int argc, char** argv){
	ElementT<branch> b;
	Problem p;
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