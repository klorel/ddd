#include "MomentGenerator.h"
#include "ComplexPolynomial.h"

MomentGenerator::MomentGenerator(size_t nvariables, size_t order):_nvariables(nvariables), _order(order) {

}

MomentGenerator::~MomentGenerator() {

}
//std::ostream & MomentGenerator::print(std::ostream & stream, Int2Int const & alpha, size_t n) const{
//	int i(0);
//	for (auto const & kvp : alpha){
//		while (i < kvp.first){
//			stream << 0;
//			++i;
//		}
//		stream << kvp.second;
//		i = kvp.first;
//		++i;
//	}
//	while (i < n){
//		stream << 0;
//		++i;
//	}
//	return stream;
//}

void MomentGenerator::build(ComplexMonomialPtrList & result) {
	result.clear();
	result.push_back(ComplexMonomial::ZeroPtr);
	int p(0);
	while (p < _nvariables) {
		ComplexMonomialPtrList new_alpha_sum;
		ComplexMonomialPtr generator(ComplexMonomial::Build(p));
		for (auto const & kvp : result) {
			int const degree(kvp->degree());
			ComplexMonomialPtr temp(kvp);
			new_alpha_sum.push_back(temp);
			for (int i(0); i < _order - degree;++i){
				temp = *temp+*generator;
				new_alpha_sum.push_back(temp);
			}
		}	
		result = new_alpha_sum;
		++p;
	}
	//std::cout << "new_alpha_sum is " << std::endl;
	//for (auto const & kvp : alphaSums)
	//	std::cout << kvp << std::endl;
	//	printAlpha(std::cout << kvp->second << " | ",kvp->first) << std::endl;
}
