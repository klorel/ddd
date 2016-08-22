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

void MomentGenerator::build() {
	std::cout << "building moment of " << _nvariables << " variables, order is " << _order << std::endl;
	Int2Int alpha;
	_alphas.clear();
	bool stop(false);

	typedef ComplexMonomial AlphaSum;
	//typedef std::pair<IntVector, int> AlphaSum;
	typedef std::shared_ptr<AlphaSum> AlphaSumPtr;
	typedef std::list<AlphaSumPtr> AlphaSumPtrList;


	AlphaSumPtr y0ptr(new ComplexMonomial);
	//*y0ptr = { IntVector(_nvariables, 0), 0 };
	//y0ptr->second = 0;
	AlphaSumPtrList alphaSums({ y0ptr });
	int p(0);
	while (p < _nvariables) {
		AlphaSumPtrList new_alpha_sum;
		for (auto const & kvp : alphaSums) {
			ComplexPolynomial temp(kvp);
			int const degree(kvp->degree());
			for (int i(0); i <= _order - degree;++i){
				temp *= ComplexPolynomial(p);
				ComplexPolynomial const & c_temp(temp);
				new_alpha_sum.push_back(c_temp.terms().begin()->first);
				std::cout << c_temp << std::endl;
			//	*new_alpha_sum.back() = temp*ComplexPolynomial::Build(p);
			//	//printAlpha(std::cout << new_alpha_sum.back()->second << " | ", new_alpha_sum.back()->first) << std::endl;
			}
		}	
		alphaSums = new_alpha_sum;
		++p;
	}
	//std::cout << "new_alpha_sum is " << std::endl;
	//for (auto const & kvp : alphaSums)
	//	printAlpha(std::cout << kvp->second << " | ",kvp->first) << std::endl;
	std::cout << "size of moment is " << alphaSums.size() << std::endl;
}
