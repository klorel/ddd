#include "MomentGenerator.h"


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


void MomentGenerator::build(){
	std::cout << "building moment of " << _nvariables << " variables, order is " << _order << std::endl;
	Int2Int alpha;
	_alphas.clear();
	bool stop(false);
	typedef std::pair<IntVector, int> AlphaSum;
	typedef std::list<AlphaSum> AlphaSums;
	AlphaSums alphaSums({ { IntVector(_nvariables, 0), 0 } });
	int p(0);
	while (p<_nvariables){		
		AlphaSums new_alpha_sum;
		for (auto const & kvp : alphaSums){
			for (int i(0); i <= _order - kvp.second;++i){
				new_alpha_sum.push_back(kvp);
				new_alpha_sum.back().first[p] = i;
				new_alpha_sum.back().second += i;
			}
		}	
		alphaSums = new_alpha_sum;
		++p;
	}
	//std::cout << "new_alpha_sum is " << std::endl;
	//for (auto const & kvp : alphaSums)
	//	printAlpha(std::cout << kvp.second << " | ", kvp.first) << std::endl;
	std::cout << "size of moment is " << alphaSums.size() << std::endl;
}
